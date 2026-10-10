og_gauss_hermite <- function(nq) {
  if (nq == 1) return(list(x = 0, w = 1))
  J <- matrix(0, nq, nq)
  off <- sqrt(1:(nq - 1))
  J[cbind(1:(nq - 1), 2:nq)] <- off
  J[cbind(2:nq, 1:(nq - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}
