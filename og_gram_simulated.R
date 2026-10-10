og_gram_simulated <- function(mom, alphas) {
  P <- og_features(mom, alphas)
  crossprod(P) / nrow(P)
}