og_simulate <- function(n, process, design, sigma2 = 0, gfun, r, oracle = TRUE, nq = 20) {
  L <- process$L
  if (r > L) stop("g may depend only on the first L latent scores")
  times <- design(n)
  V <- matrix(rnorm(n * L), n, L)
  w <- vector("list", n)
  for (i in seq_len(n)) {
    tt <- times[[i]]; ri <- length(tt)
    if (ri == 0) { w[[i]] <- numeric(0); next }
    A <- process$cross(tt)
    Ct <- process$cov(tt, tt)
    x <- drop(A %*% V[i, ]) + drop(og_psd_sqrt(Ct - A %*% t(A), scale = max(abs(Ct))) %*% rnorm(ri))
    w[[i]] <- if (is.function(sigma2)) x + drop(og_psd_sqrt(sigma2(tt)) %*% rnorm(ri)) else
      x + sqrt(sigma2) * rnorm(ri)
  }
  g <- gfun(V[, seq_len(r), drop = FALSE])
  if (any(!is.finite(g)) || any(g < 0 | g > 1)) stop("gfun must return probabilities in [0, 1]")
  Y <- rbinom(n, 1, g)
  out <- list(t = times, w = w, Y = Y, V = V, g = g)
  if (oracle) out$Kg <- og_oracle_Kg(og_moments_all(out, process, sigma2), gfun, r, nq)
  out
}
