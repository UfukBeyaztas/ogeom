og_bm_augmented <- function(T, W, Y, J, j_all = og_bessel_j0_zeros(J + 2000)) {
  if (any(T < 0 | T > 1)) stop("times must lie in [0, 1]")
  if (any(T == 0 & W != 0)) stop("a measurement at time 0 must equal 0 for error-free Brownian motion")
  n <- length(Y)
  tails <- og_fb_tails(J, j_all)
  rJ <- sqrt(tails[["t1"]]); mu <- tails[["t2"]] / tails[["t1"]]
  j <- j_all[seq_len(J)]
  one_k <- 2 * sign(besselJ(j, 1)) / j          # <1, phi_k>
  pos <- T > 0
  Psi <- matrix(0, n, J + 1)
  if (any(pos)) {
    Tp <- T[pos]; Wp <- W[pos]
    cumk <- if (J > 0) og_fb_cum(Tp, j) else matrix(0, length(Tp), 0)
    proj <- if (J > 0) drop(cumk %*% one_k) else rep(0, length(Tp))
    if (J > 0) Psi[pos, seq_len(J)] <- cumk / Tp * Wp
    Psi[pos, J + 1] <- Wp / (Tp * rJ) * (Tp - proj)   # T^{-1} <e_J, 1_[0,T]> W
  }
  chat <- colMeans((Y - 0.5) * Psi)
  theta <- chat / c(4 / j^2, mu)
  list(theta = theta[seq_len(J)], theta_e = theta[J + 1], j = j, lambda = 4 / j^2,
       mu = mu, rJ = rJ, one_k = one_k)
}
