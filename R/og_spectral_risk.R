og_spectral_risk <- function(lambda, theta, n, tail) {
  sum(1 / (4 * n * lambda) - theta^2 / n) + tail
}
