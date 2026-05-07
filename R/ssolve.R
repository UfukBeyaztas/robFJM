ssolve <- function(A, b, ridge = 1e-8) {
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(out))) return(out)
  out2 <- tryCatch(solve(A + diag(ridge, 2), b), error = function(e) c(0,0))
  out2[!is.finite(out2)] <- 0
  out2
}
