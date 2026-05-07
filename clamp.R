clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}