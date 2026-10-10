og_psd_sqrt <- function(S, scale = NULL, tol = 1e-8) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  if (is.null(scale)) scale <- max(abs(e$values))
  if (length(e$values) > 0 && min(e$values) < -tol * scale)
    stop("og_psd_sqrt: matrix is not positive semidefinite at the chosen tolerance")
  e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}
