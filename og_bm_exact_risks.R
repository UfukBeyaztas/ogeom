og_bm_exact_risks <- function(c0, theta0, J, n, j_all) {
  K <- length(j_all)
  if (length(theta0) > K) stop("theta0 is longer than j_all")
  lam <- 4 / j_all^2; one <- 2 * sign(besselJ(j_all, 1)) / j_all
  th0 <- numeric(K); th0[seq_along(theta0)] <- theta0
  th <- c0 * one + th0
  tl <- og_fb_tails(J, j_all)
  idx <- if (J < K) (J + 1):K else integer(0)
  head <- seq_len(J)
  tail_sq <- c0^2 * tl[["t1"]] + 2 * c0 * sum(one[idx] * th0[idx]) + sum(th0[idx]^2)
  spectral <- sum(1 / (4 * n * lam[head]) - th[head]^2 / n) + tail_sq
  rJ <- sqrt(tl[["t1"]]); mu <- tl[["t2"]] / tl[["t1"]]
  theta_e <- (c0 * tl[["t2"]] + sum(lam[idx] * one[idx] * th0[idx])) / (rJ * mu)
  t_e <- (c0 * tl[["t1"]] + sum(one[idx] * th0[idx])) / rJ
  bias <- tail_sq - 2 * theta_e * t_e + theta_e^2
  augmented <- (sum(1 / lam[head]) + 1 / mu) / (4 * n) - (sum(th[head]^2) + theta_e^2) / n + bias
  c(spectral = spectral, augmented = augmented, tail_sq = tail_sq, theta_e = theta_e, t_e = t_e)
}