og_fb_cum <- function(t, j) {
  X <- outer(sqrt(t), j)
  M <- 2 * sqrt(t) * matrix(besselJ(as.vector(X), 1), nrow(X))
  sweep(M, 2, j * abs(besselJ(j, 1)), "/")
}