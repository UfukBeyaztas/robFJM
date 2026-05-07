init_params <- function(dat_surv, dat_long, long_formula, surv_formula) {
  lm_fit <- lm(long_formula, data = dat_long)
  beta0 <- coef(lm_fit); beta0[is.na(beta0)] <- 0
  sig_e2 <- max(mean(resid(lm_fit)^2), 1e-3)
  
  sr <- survreg(surv_formula, dist = "weibull", data = dat_surv)
  b_aft <- coef(sr); sc <- sr$scale
  Tau0 <- 1/sc
  gamma0 <- -b_aft/sc
  
  list(beta = beta0, gamma = gamma0, Tau = Tau0, phi = 0.1,
       sigmasq.e = sig_e2, sigmasq.u = 0.5)
}