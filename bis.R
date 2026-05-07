bis <- function(x, tol) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(FALSE)
  isTRUE(max(abs(x)) < tol)
}