og_oracle_Kg <- function(mom, gfun, r, nq = 20) {
  gh <- og_gauss_hermite(nq)
  Zg <- as.matrix(expand.grid(rep(list(gh$x), r)))
  wts <- apply(as.matrix(expand.grid(rep(list(gh$w), r))), 1, prod)
  n <- nrow(mom$Uhat); out <- numeric(n)
  for (i in seq_len(n)) {
    S <- diag(r) - matrix(mom$beta[i, 1:r, 1:r], r, r)
    pts <- sweep(Zg %*% og_psd_sqrt(S, scale = 1), 2, mom$Uhat[i, 1:r], "+")
    out[i] <- sum(wts * gfun(pts))
  }
  out
}