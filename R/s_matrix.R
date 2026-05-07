s_matrix <- function(M, intercept_name = "(Intercept)", protect = NULL) {
  cols <- colnames(M)
  is_int <- cols == intercept_name
  is_protect <- rep(FALSE, ncol(M))
  if (!is.null(protect)) is_protect <- cols %in% protect
  keep_raw <- is_int | is_protect
  
  mu <- rep(0, ncol(M))
  sd <- rep(1, ncol(M))
  
  if (any(!keep_raw)) {
    mu[!keep_raw] <- colMeans(M[, !keep_raw, drop = FALSE])
    sd[!keep_raw] <- apply(M[, !keep_raw, drop = FALSE], 2, sd)
    sd[!keep_raw] <- ifelse(!is.finite(sd[!keep_raw]) | sd[!keep_raw] == 0, 1, sd[!keep_raw])
    M[, !keep_raw] <- sweep(M[, !keep_raw, drop = FALSE], 2, mu[!keep_raw], "-")
    M[, !keep_raw] <- sweep(M[, !keep_raw, drop = FALSE], 2, sd[!keep_raw], "/")
  }
  
  list(M = M, mu = mu, sd = sd, is_int = is_int, cols = cols, keep_raw = keep_raw)
}
