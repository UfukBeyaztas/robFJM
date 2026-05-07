fit_APHL <- function(des, init,
                     nu = 4,
                     max_outer = 30, inner_u = 5, inner_th = 10,
                     tol = 1e-6,
                     Tau_bounds = c(0.2, 3.0),
                     phi_bounds = c(-2, 2),
                     sig_bounds = c(1e-6, 1e6),
                     u_bounds = c(-8, 8),
                     damp_u = 0.5, damp_pa = 0.25, damp_beta = 1.0,
                     exp_cap = 25) {
  
  Xcols <- des$X_scale$cols
  Zcols <- des$Z_scale$cols
  time_col <- "obstime"
  
  X <- des$X_all
  Z <- des$Z
  
  time <- as.numeric(des$time)
  d <- as.numeric(des$event)
  
  n_long <- nrow(des$dat_long)
  k <- length(des$ids)
  p <- ncol(X)
  
  time_idx <- match(time_col, colnames(X))
  if (is.na(time_idx)) stop("fit_APHL: cannot find 'obstime' column in X.")
  
  beta <- init$beta[Xcols]; beta[is.na(beta)] <- 0
  gamma <- init$gamma[Zcols]; gamma[is.na(gamma)] <- 0
  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)
  
  sig_e2 <- clamp(init$sigmasq.e, sig_bounds[1], sig_bounds[2])
  sig_u2 <- clamp(init$sigmasq.u, sig_bounds[1], sig_bounds[2])
  
  a <- log(clamp(init$Tau, Tau_bounds[1], Tau_bounds[2]))
  phi <- clamp(init$phi, phi_bounds[1], phi_bounds[2])
  
  u <- rep(0, k)
  w_all <- rep(1, n_long)
  
  compute_quadrature <- function(b_vec, g_vec, u_vec, a_val, phi_val) {
    Tau_val <- exp(a_val)
    eta_S_val <- drop(Z %*% g_vec)
    
    out <- list(
      Int_g = numeric(k), Int_g_m = numeric(k), Int_g_tau = numeric(k),
      Int_g_m2 = numeric(k), Int_g_aa = numeric(k), Int_g_m_tau = numeric(k),
      m_Ti = numeric(k),
      Int_g_D = matrix(0, nrow = k, ncol = p),
      Int_g_DD = array(0, dim = c(p, p, k)),
      D_Ti = matrix(0, nrow = k, ncol = p)
    )
    
    for (i in 1:k) {
      T_i <- max(time[i], 1e-6)
      Xi_base <- des$Xlist[[i]][1, ]
      Xi_base <- as.numeric(Xi_base)
      
      N_grid <- 30
      t_grid <- seq(1e-6, T_i, length.out = N_grid)
      
      dx <- diff(t_grid)
      wts <- numeric(N_grid)
      wts[1] <- dx[1] / 2
      wts[N_grid] <- dx[N_grid - 1] / 2
      if (N_grid > 2) wts[2:(N_grid - 1)] <- (dx[1:(N_grid - 2)] + dx[2:(N_grid - 1)]) / 2
      
      D_grid <- matrix(rep(Xi_base, each = N_grid), nrow = N_grid)
      D_grid[, time_idx] <- t_grid
      
      m_grid <- drop(D_grid %*% b_vec) + u_vec[i]
      out$m_Ti[i] <- m_grid[N_grid]
      out$D_Ti[i, ] <- D_grid[N_grid, ]
      
      g_grid <- Tau_val * (t_grid^(Tau_val - 1)) *
        sexp(eta_S_val[i] + phi_val * m_grid, lim = exp_cap)
      
      tau_logt <- Tau_val * log(t_grid)
      term_tau <- 1 + tau_logt
      term_aa  <- term_tau^2 + tau_logt
      
      wg <- wts * g_grid
      
      out$Int_g[i]       <- sum(wg)
      out$Int_g_m[i]     <- sum(wg * m_grid)
      out$Int_g_tau[i]   <- sum(wg * term_tau)
      out$Int_g_m2[i]    <- sum(wg * m_grid^2)
      out$Int_g_aa[i]    <- sum(wg * term_aa)
      out$Int_g_m_tau[i] <- sum(wg * m_grid * term_tau)
      
      out$Int_g_D[i, ]   <- colSums(D_grid * wg)
      
      weighted_D <- D_grid * sqrt(pmax(wg, 0))
      out$Int_g_DD[,,i]  <- crossprod(weighted_D)
    }
    
    out
  }
  
  # ----------------- main loop -----------------
  for (outer in 1:max_outer) {
    par_old <- c(beta, gamma, a, phi, sig_e2, sig_u2, u)
    
    # ---- 1) E-step weights
    u_long <- u[match(des$dat_long$id, des$ids)]
    r_all <- des$dat_long$yij - drop(X %*% beta) - u_long
    w_all <- (nu + 1) / (nu + (r_all^2) / pmax(sig_e2, 1e-12))
    w_all[!is.finite(w_all)] <- 0
    w_all <- pmax(w_all, 0)
    
    quad <- compute_quadrature(beta, gamma, u, a, phi)
    
    # ---- 2) update u (Newton)
    for (it in 1:inner_u) {
      score_u <- numeric(k)
      info_u  <- numeric(k)
      
      for (i in 1:k) {
        idx <- des$idx_by_id[[i]]
        wi <- w_all[idx]
        ri <- des$ylist[[i]] - drop(des$Xlist[[i]] %*% beta) - u[i]
        
        long_score <- sum(wi * ri) / pmax(sig_e2, 1e-12)
        long_info  <- sum(wi) / pmax(sig_e2, 1e-12)
        
        score_u[i] <- long_score + phi * (d[i] - quad$Int_g[i]) - u[i] / pmax(sig_u2, 1e-12)
        info_u[i]  <- long_info  + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
      }
      
      step <- score_u / pmax(info_u, 1e-12)
      step[!is.finite(step)] <- 0
      u <- clamp(u + damp_u * step, u_bounds[1], u_bounds[2])
      
      quad <- compute_quadrature(beta, gamma, u, a, phi)
      if (bis(step, 1e-8)) break
    }
    
    # ---- 3) refresh residuals/weights and Hinv
    u_long <- u[match(des$dat_long$id, des$ids)]
    r_all <- des$dat_long$yij - drop(X %*% beta) - u_long
    w_all <- (nu + 1) / (nu + (r_all^2) / pmax(sig_e2, 1e-12))
    w_all[!is.finite(w_all)] <- 0
    w_all <- pmax(w_all, 0)
    
    H_u <- numeric(k)
    for (i in 1:k) {
      idx <- des$idx_by_id[[i]]
      H_u[i] <- (sum(w_all[idx]) / pmax(sig_e2, 1e-12)) + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
    }
    Hinv <- 1 / pmax(H_u, 1e-12)
    
    # ---- 4) beta update (APHL Newton)
    U_beta_L <- drop(crossprod(X, w_all * r_all)) / pmax(sig_e2, 1e-12)
    
    U_beta_S <- rep(0, p)
    U_beta_trace <- rep(0, p)
    I_beta_S <- matrix(0, p, p)
    
    for (i in 1:k) {
      U_beta_S <- U_beta_S + d[i] * phi * quad$D_Ti[i, ] - phi * quad$Int_g_D[i, ]
      U_beta_trace <- U_beta_trace - 0.5 * Hinv[i] * (phi^3) * quad$Int_g_D[i, ]
      I_beta_S <- I_beta_S + (phi^2) * quad$Int_g_DD[,,i]
    }
    
    U_beta <- U_beta_L + U_beta_S + U_beta_trace
    I_beta_L <- crossprod(X * sqrt(w_all)) / pmax(sig_e2, 1e-12)
    I_beta <- I_beta_L + I_beta_S + diag(1e-8, p)
    
    step_beta <- tryCatch(solve(I_beta, U_beta), error = function(e) rep(0, p))
    step_beta[!is.finite(step_beta)] <- 0
    beta <- beta + damp_beta * step_beta
    names(beta) <- colnames(X)
    
    u_long <- u[match(des$dat_long$id, des$ids)]
    r_all <- des$dat_long$yij - drop(X %*% beta) - u_long
    w_all <- (nu + 1) / (nu + (r_all^2) / pmax(sig_e2, 1e-12))
    w_all[!is.finite(w_all)] <- 0
    w_all <- pmax(w_all, 0)
    
    quad <- compute_quadrature(beta, gamma, u, a, phi)
    for (i in 1:k) {
      idx <- des$idx_by_id[[i]]
      H_u[i] <- (sum(w_all[idx]) / pmax(sig_e2, 1e-12)) + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
    }
    Hinv <- 1 / pmax(H_u, 1e-12)
    
    # ---- 5) update gamma, phi, logTau (APHL)
    for (it in 1:inner_th) {
      
      # gamma step
      r_surv_adj <- d - quad$Int_g - 0.5 * Hinv * (phi^2) * quad$Int_g
      S_g <- drop(crossprod(Z, r_surv_adj))
      Hdiag_g <- pmax(drop(crossprod(Z^2, quad$Int_g)), 1e-12)
      step_g <- S_g / Hdiag_g
      step_g[!is.finite(step_g)] <- 0
      gamma_old <- gamma
      gamma <- gamma + 0.5 * step_g
      names(gamma) <- colnames(Z)
      
      quad <- compute_quadrature(beta, gamma, u, a, phi)
      for (i in 1:k) {
        idx <- des$idx_by_id[[i]]
        H_u[i] <- (sum(w_all[idx]) / pmax(sig_e2, 1e-12)) + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
      }
      Hinv <- 1 / pmax(H_u, 1e-12)
      
      Tau <- exp(a)
      logT <- log(pmax(time, 1e-6))
      
      S_phi <- sum(d * quad$m_Ti - quad$Int_g_m -
                     0.5 * Hinv * (2 * phi * quad$Int_g + (phi^2) * quad$Int_g_m))
      
      S_a <- sum(d * (1 + Tau * logT) - quad$Int_g_tau -
                   0.5 * Hinv * (phi^2) * quad$Int_g_tau)
      
      I_pp <- sum(quad$Int_g_m2) + 1e-8
      I_aa <- sum(-d * Tau * logT + quad$Int_g_aa) + 1e-8
      I_pa <- sum(quad$Int_g_m_tau)
      
      step2 <- ssolve(matrix(c(I_pp, I_pa, I_pa, I_aa), 2, 2), c(S_phi, S_a))
      step2[!is.finite(step2)] <- 0
      
      phi <- clamp(phi + damp_pa * step2[1], phi_bounds[1], phi_bounds[2])
      a   <- clamp(a   + damp_pa * step2[2], log(Tau_bounds[1]), log(Tau_bounds[2]))
      
      quad <- compute_quadrature(beta, gamma, u, a, phi)
      for (i in 1:k) {
        idx <- des$idx_by_id[[i]]
        H_u[i] <- (sum(w_all[idx]) / pmax(sig_e2, 1e-12)) + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
      }
      Hinv <- 1 / pmax(H_u, 1e-12)
      
      if (bis(c(gamma - gamma_old, step2), 1e-6)) break
    }
    
    # ---- 6) variances 
    sig_e2_new <- sum(w_all * r_all^2) / max(n_long, 1)
    sig_e2 <- clamp(sig_e2_new, sig_bounds[1], sig_bounds[2])
    
    for (i in 1:k) {
      idx <- des$idx_by_id[[i]]
      H_u[i] <- (sum(w_all[idx]) / pmax(sig_e2, 1e-12)) + (phi^2) * quad$Int_g[i] + 1 / pmax(sig_u2, 1e-12)
    }
    H_u <- pmax(H_u, 1e-12)
    
    sig_u2_new <- mean(u^2 + 1 / H_u)
    sig_u2 <- clamp(sig_u2_new, sig_bounds[1], sig_bounds[2])
    
    par_new <- c(beta, gamma, a, phi, sig_e2, sig_u2, u)
    if (bis(par_new - par_old, tol)) break
  }
  
  beta_u  <- uns_coef(beta,  des$X_scale$mu, des$X_scale$sd, des$X_scale$is_int, des$X_scale$keep_raw)
  gamma_u <- uns_coef(gamma, des$Z_scale$mu, des$Z_scale$sd, des$Z_scale$is_int, des$Z_scale$keep_raw)
  
  list(beta = beta_u, gamma = gamma_u, Tau = exp(a), phi = phi,
       sig_e2 = sig_e2, sig_u2 = sig_u2, u_hat = u)
}