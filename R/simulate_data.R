simulate_data <- function(I = 200, J = 4,
                          D = 10, by = 0.1,
                          beta0 = 1.35, beta1 = 0.78,
                          var_b = 1, var_eps = 1,
                          gamma0 = 0, gamma1 = -1.75, alpha = 0.29,
                          Tau = 0.8,
                          lambdaC = -6.0,
                          Kx_latent = 5, Kw_latent = 5,
                          meas_sd_X = 0.15, meas_sd_W = 0.15,
                          outlier_frac_Y = 0.00, outlier_mag_Y = c(20, 30),
                          outlier_frac_fun = 0.00, outlier_mag_fun = c(3, 6)) {
  
  s <- seq(0, D, by)
  m <- length(s)
  
  trueBeta  <- 2.0 * sin(s * pi / 5)
  trueGamma <- 1.2 * sin(s * pi / 4)
  
  Phi_latent <- sapply(1:Kx_latent, function(k) sin(k * pi * s / D))
  Psi_latent <- sapply(1:Kw_latent, function(k) cos(k * pi * s / D))
  
  norm_col <- function(A) {
    for (k in 1:ncol(A)) {
      nk <- sqrt(trapz(s, A[,k]^2))
      if (!is.finite(nk) || nk == 0) nk <- 1
      A[,k] <- A[,k] / nk
    }
    A
  }
  Phi_latent <- norm_col(Phi_latent)
  Psi_latent <- norm_col(Psi_latent)
  
  W1 <- rbinom(I, 1, 0.5)
  b_i <- rnorm(I, 0, sqrt(var_b))
  
  xi <- matrix(rnorm(I*Kx_latent), I, Kx_latent) %*% diag(1/(1:Kx_latent))
  zi <- matrix(rnorm(I*Kw_latent), I, Kw_latent) %*% diag(1/(1:Kw_latent))
  
  # Latent curves
  X_true <- xi %*% t(Phi_latent)
  W_true <- zi %*% t(Psi_latent)
  
  # Observed curves
  Xmat <- X_true + matrix(rnorm(I*m, 0, meas_sd_X), I, m)
  Wmat <- W_true + matrix(rnorm(I*m, 0, meas_sd_W), I, m)
  
  # ---- FUNCTIONAL OUTLIERS (Magnitude Shifts) ----
  if (outlier_frac_fun > 0) {
    n_out <- max(1, floor(outlier_frac_fun * I))
    id_out <- sample.int(I, n_out)
    
    for (ii in id_out) {
      shift_X <- runif(1, outlier_mag_fun[1], outlier_mag_fun[2]) * sample(c(-1,1), 1)
      shift_W <- runif(1, outlier_mag_fun[1], outlier_mag_fun[2]) * sample(c(-1,1), 1)
      
      Xmat[ii, ] <- Xmat[ii, ] + shift_X
      Wmat[ii, ] <- Wmat[ii, ] + shift_W
    }
  }
  
  intX <- apply(X_true, 1, function(x) trapz(s, x * trueBeta))
  intW <- apply(W_true, 1, function(w) trapz(s, w * trueGamma))
  
  obstime_vals <- seq(0, D, length.out = J)
  
  eta_surv_base <- gamma0 + gamma1*W1 + intW
  U <- runif(I)
  Ttrue <- numeric(I)
  
  for (i in 1:I) {
    hazard_func <- function(t) {
      m_t <- beta0 + beta1 * t + intX[i] + b_i[i]
      Tau * (t^(Tau - 1)) * exp(eta_surv_base[i] + alpha * m_t)
    }
    
    target_func <- function(Time) {
      if (Time <= 1e-6) return(0 - (-log(U[i])))
      int_val <- tryCatch(integrate(hazard_func, lower = 1e-6, upper = Time)$value, error = function(e) NA)
      if (is.na(int_val)) return(NA)
      return(int_val + log(U[i]))
    }
    
    root_res <- tryCatch(
      uniroot(target_func, interval = c(1e-5, 200), extendInt = "yes"),
      error = function(e) list(root = 200) 
    )
    Ttrue[i] <- root_res$root
  }
  
  Uc <- runif(I)
  C <- (-log(Uc) / exp(lambdaC))^(1/Tau)
  C <- pmin(C, 180)
  
  time <- pmin(Ttrue, C)
  event <- as.integer(Ttrue <= C)
  
  surv <- data.frame(ID = 1:I, time = time, event = event, W1 = W1)
  colnames(Wmat) <- paste0("func.W.", 1:m)
  surv <- cbind(surv, as.data.frame(Wmat))
  
  long_list <- vector("list", I)
  colnames(Xmat) <- paste0("func.X.", 1:m)
  
  for (i in 1:I) {
    eps_ij <- rnorm(J, 0, sqrt(var_eps))
    Y <- beta0 + beta1*obstime_vals + intX[i] + b_i[i] + eps_ij
    
    df_i <- data.frame(ID = i, visit = 1:J, obstime = obstime_vals, time = time[i], event = event[i], Y = Y)
    df_i <- cbind(df_i, as.data.frame(matrix(rep(Xmat[i,], each = J), nrow = J, byrow = FALSE)))
    colnames(df_i)[(ncol(df_i)-m+1):ncol(df_i)] <- paste0("func.X.", 1:m)
    
    df_i$b <- b_i[i]
    long_list[[i]] <- df_i
  }
  long <- do.call(rbind, long_list)
  
  # ---- LONGITUDINAL OUTLIERS (Symmetric Heavy-Tailed) ----
  if (outlier_frac_Y > 0) {
    n_obs <- nrow(long)
    m_out <- max(1, floor(outlier_frac_Y * n_obs))
    idx <- sample.int(n_obs, m_out)
    long$Y[idx] <- long$Y[idx] + runif(m_out, outlier_mag_Y[1], outlier_mag_Y[2]) * sample(c(-1,1), m_out, TRUE)
  }
  
  trimmed <- tltf(long, surv)
  
  list(
    data.long = trimmed$long,
    data.surv = trimmed$surv,
    grid = s,
    by = by,
    true = list(
      beta0 = beta0, beta1 = beta1,
      var_b = var_b, var_eps = var_eps,
      gamma0 = gamma0, gamma1 = gamma1,
      alpha = alpha, Tau = Tau,
      trueBeta = trueBeta, trueGamma = trueGamma
    )
  )
}
