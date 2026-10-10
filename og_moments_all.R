og_moments_all <- function(data, process, sigma2 = 0, tol = 1e-10) {
  n <- length(data$t); L <- process$L
  Uhat <- matrix(0, n, L); beta <- array(0, c(n, L, L))
  for (i in seq_len(n)) {
    m <- og_conditional_moments(data$t[[i]], data$w[[i]], process, sigma2, tol)
    Uhat[i, ] <- m$Uhat
    beta[i, , ] <- m$beta
  }
  list(Uhat = Uhat, beta = beta)
}