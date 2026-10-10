og_bm_spectral <- function(T, W, Y, J, j = og_bessel_j0_zeros(J)) {
  if (any(T < 0 | T > 1)) stop("times must lie in [0, 1]")
  if (any(T == 0 & W != 0)) stop("a measurement at time 0 must equal 0 for error-free Brownian motion")
  if (J == 0) return(list(theta = numeric(0), j = numeric(0), lambda = numeric(0)))
  j <- j[seq_len(J)]
  pos <- T > 0
  Psi <- matrix(0, length(T), J)
  if (any(pos)) Psi[pos, ] <- og_fb_cum(T[pos], j) / T[pos] * W[pos]
  chat <- colMeans((Y - 0.5) * Psi)
  list(theta = chat * j^2 / 4, j = j, lambda = 4 / j^2)
}
