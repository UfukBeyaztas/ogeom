og_bm_linear_risks <- function(c0, c1, J, n, j_all) {
  if (length(J) != 1 || J < 0 || J != as.integer(J) || J > min(100, length(j_all)))
    stop("J must be an integer between 0 and min(100, length(j_all))")
  if (n < 1) stop("n must be positive")
  head <- seq_len(J)
  lam <- 4 / j_all[head]^2
  one <- 2 * sign(besselJ(j_all[head], 1)) / j_all[head]
  th <- (c0 + c1 * lam) * one
  tl <- c(t1 = 1 - sum(lam), t2 = 0.5 - sum(lam^2), t3 = 1/3 - sum(lam^3))
  tail_sq <- c0^2 * tl[["t1"]] + 2 * c0 * c1 * tl[["t2"]] + c1^2 * tl[["t3"]]
  r <- sqrt(tl[["t1"]]); mu <- tl[["t2"]] / tl[["t1"]]
  t_e <- (c0 * tl[["t1"]] + c1 * tl[["t2"]]) / r
  theta_e <- (c0 * tl[["t2"]] + c1 * tl[["t3"]]) / (r * mu)
  bias <- max(0, tail_sq - 2 * theta_e * t_e + theta_e^2)
  c(spectral = sum(1 / (4 * n * lam) - th^2 / n) + tail_sq,
    augmented = (sum(1 / lam) + 1 / mu) / (4 * n) - (sum(th^2) + theta_e^2) / n + bias,
    tail_sq = tail_sq, theta_e = theta_e, t_e = t_e)
}