uns_coef <- function(b_std, mu, sd, is_int, keep_raw) {
  b <- b_std
  idx_std <- !keep_raw
  if (any(idx_std)) {
    b[idx_std] <- b_std[idx_std] / sd[idx_std]
    b0 <- b_std[is_int]
    adj <- sum((b_std[idx_std] * mu[idx_std]) / sd[idx_std])
    b[is_int] <- b0 - adj
  }
  b
}