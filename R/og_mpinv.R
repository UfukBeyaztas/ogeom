og_mpinv <- function(G, tol = 1e-10) {
  G <- (G + t(G)) / 2
  e <- eigen(G, symmetric = TRUE)
  scale <- max(abs(e$values))
  if (scale == 0) return(matrix(0, nrow(G), ncol(G)))
  if (min(e$values) < -tol * scale) stop("og_mpinv: matrix is not positive semidefinite")
  keep <- e$values > tol * scale
  V <- e$vectors[, keep, drop = FALSE]
  V %*% (t(V) / e$values[keep])
}
