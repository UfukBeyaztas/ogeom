og_ure <- function(eta, lambda, n) {
  if (n < 2) stop("the risk criterion needs n >= 2")
  c(0, cumsum(1 / (2 * (n - 1) * lambda) - (n + 1) / (n - 1) * eta^2))
}
