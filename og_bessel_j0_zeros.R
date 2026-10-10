og_bessel_j0_zeros <- function(K) {
  if (K == 0) return(numeric(0))
  vapply(seq_len(K), function(k)
    uniroot(function(x) besselJ(x, 0), c((k - 0.25) * pi, k * pi), tol = 1e-13)$root,
    numeric(1))
}