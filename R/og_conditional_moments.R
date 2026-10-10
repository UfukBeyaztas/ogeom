og_conditional_moments <- function(times, w, process, sigma2 = 0, tol = 1e-10) {
  L <- process$L; r <- length(times)
  if (r == 0) return(list(Uhat = numeric(L), beta = matrix(0, L, L)))
  A <- process$cross(times)
  Sig <- if (is.function(sigma2)) sigma2(times) else diag(sigma2, r)
  G <- process$cov(times, times) + Sig
  GA <- og_mpinv(G, tol) %*% A
  beta <- crossprod(A, GA)
  list(Uhat = drop(crossprod(GA, w)), beta = (beta + t(beta)) / 2)
}
