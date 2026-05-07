fit_FAPHL <- function(long, surv, FPCA.X, FPCA.W, by,
                      nu = 4, verbose = FALSE) {
  
  dat_long <- long
  dat_surv <- surv
  
  dat_long$id <- dat_long$ID
  dat_long$yij <- dat_long$Y
  
  dat_surv$id <- dat_surv$ID
  dat_surv$si <- dat_surv$time
  
  long_formula <- yij ~ obstime + FPCscore.X
  surv_formula <- Surv(si, event) ~ W1 + FPCscore.W
  
  des <- build_designs_std(dat_surv, dat_long, long_formula, surv_formula)
  init <- init_params(des$dat_surv, des$dat_long, long_formula, surv_formula)
  
  fit <- fit_APHL(des, init, nu = nu, max_outer = 30)
  
  bX <- extract_pc_coefs(fit$beta,  prefix = "FPCscore.X")
  gW <- extract_pc_coefs(fit$gamma, prefix = "FPCscore.W")
  
  BetaHat  <- if (length(bX) == ncol(FPCA.X$basis)) as.numeric(FPCA.X$basis %*% bX) else rep(NA_real_, nrow(FPCA.X$basis))
  GammaHat <- if (length(gW) == ncol(FPCA.W$basis)) as.numeric(FPCA.W$basis %*% gW) else rep(NA_real_, nrow(FPCA.W$basis))
  
  scalars <- c(
    beta1   = unname(fit$beta["obstime"]),
    var_b   = unname(fit$sig_u2),
    var_eps = unname(fit$sig_e2),
    gamma1  = unname(fit$gamma["W1"]),
    alpha   = unname(fit$phi),
    Tau     = unname(fit$Tau)
  )
  
  list(fit = fit, BetaHat = BetaHat, GammaHat = GammaHat, scalars = scalars)
}
