og_risk_exact <- function(A, G, c, n, theta_target) {
  m <- drop(A %*% c)
  sum(diag(A %*% G %*% t(A))) / (4 * n) - sum(m^2) / n + sum((m - theta_target)^2)
}