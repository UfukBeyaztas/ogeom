og_true_coefficients <- function(gfun, r, alphas, nq = 30) {
  alphas <- as.matrix(alphas)
  gh <- og_gauss_hermite(nq)
  grid <- as.matrix(expand.grid(rep(list(gh$x), r)))
  wts <- apply(as.matrix(expand.grid(rep(list(gh$w), r))), 1, prod)
  gv <- gfun(grid)
  out <- numeric(nrow(alphas))
  for (k in seq_len(nrow(alphas))) {
    if (ncol(alphas) > r && any(alphas[k, -(1:r)] > 0)) next
    out[k] <- sum(wts * gv * og_hermite_features(grid, alphas[k, 1:r, drop = FALSE]))
  }
  names(out) <- og_alpha_labels(alphas)
  out
}
