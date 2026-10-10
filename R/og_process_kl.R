og_process_kl <- function(mu, phi, U = diag(length(mu))) {
  K <- length(mu); U <- as.matrix(U)
  if (any(!is.finite(mu)) || any(mu <= 0)) stop("eigenvalues mu must be positive and finite")
  if (nrow(U) != K) stop("U must have length(mu) rows")
  if (max(abs(crossprod(U) - diag(ncol(U)))) > 1e-8) stop("columns of U must be orthonormal")
  Phi <- function(t) matrix(phi(t), nrow = length(t), ncol = K)
  list(type = "kl", K = K, L = ncol(U), mu = mu, phi = phi, U = U,
       cov   = function(s, t) Phi(s) %*% (mu * t(Phi(t))),
       cross = function(t) Phi(t) %*% (sqrt(mu) * U))
}
