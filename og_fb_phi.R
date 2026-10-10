og_fb_phi <- function(s, j) {
  X <- outer(sqrt(s), j)
  sweep(matrix(besselJ(as.vector(X), 0), nrow(X)), 2, abs(besselJ(j, 1)), "/")
}