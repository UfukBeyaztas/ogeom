og_predict <- function(fit, newdata, clip = TRUE) {
  mom <- og_moments_all(newdata, fit$process, fit$sigma2)
  p <- rep(fit$Ybar, length(newdata$t))
  for (q in seq_len(fit$Q)) p <- p + drop(og_features(mom, fit$orders[[q]]$alphas) %*% fit$orders[[q]]$theta)
  if (clip) p <- pmin(pmax(p, 0), 1)
  p
}
