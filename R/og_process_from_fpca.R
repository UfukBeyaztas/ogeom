og_process_from_fpca <- function(fp, K = fp$selectK, sigma2_min = 0) {
  wg  <- fp$workGrid
  sgn <- sign(colMeans(fp$phi[, seq_len(K), drop = FALSE])); sgn[sgn == 0] <- 1
  phi_grid <- sweep(fp$phi[, seq_len(K), drop = FALSE], 2, sgn, "*")
  lambda   <- fp$lambda[seq_len(K)]
  phi_fun  <- function(t) vapply(seq_len(K), function(k) approx(wg, phi_grid[, k], xout = t, rule = 2)$y,
                                 numeric(length(t)))
  mu_fun   <- function(t) approx(wg, fp$mu, xout = t, rule = 2)$y
  list(process = og_process_kl(lambda, phi_fun), mu_fun = mu_fun,
       sigma2 = max(fp$sigma2, sigma2_min), lambda = lambda, phi_fun = phi_fun, sign = sgn)
}
