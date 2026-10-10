og_galerkin <- function(Psi, Y, G, J = NULL, tol = 1e-10) {
  n <- length(Y)
  e <- eigen((G + t(G)) / 2, symmetric = TRUE)
  npos <- sum(e$values > tol * max(1, e$values[1]))
  if (is.null(J)) J <- npos
  J <- min(J, npos)
  Z <- (Y - 0.5) * Psi
  chat <- colMeans(Z)
  Wj <- e$vectors[, seq_len(J), drop = FALSE]
  lam <- e$values[seq_len(J)]
  eta <- drop(crossprod(Wj, chat)) / lam
  A <- if (J > 0) Wj %*% (t(Wj) / lam) else matrix(0, ncol(Psi), ncol(Psi))
  theta <- drop(A %*% chat)
  names(theta) <- colnames(Psi)
  list(theta = theta, eta = eta, A = A, chat = chat, J = J, npos = npos,
       values = e$values, vectors = e$vectors,
       var_proxy = if (J > 0) sum(1 / lam) / (4 * n) else 0,
       cov_sandwich = A %*% cov(Z) %*% A / n)
}