og_oracle_index <- function(mom, b, link = plogis, nq = 40) {
  gh <- og_gauss_hermite(nq); r <- length(b)
  m <- drop(mom$Uhat[, 1:r, drop = FALSE] %*% b)
  bb <- 0
  for (l in 1:r) for (l2 in 1:r) bb <- bb + b[l] * b[l2] * mom$beta[, l, l2]
  s <- sqrt(pmax(sum(b^2) - bb, 0))
  drop(link(outer(m, rep(1, nq)) + outer(s, gh$x)) %*% gh$w)
}