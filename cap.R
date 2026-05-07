cap <- function(x, lim = 25) {
  pmin(pmax(x, -lim), lim)
}