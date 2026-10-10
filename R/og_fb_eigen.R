og_fb_eigen <- function(K) {
  j <- og_bessel_j0_zeros(K)
  list(j = j, lambda = 4 / j^2)
}
