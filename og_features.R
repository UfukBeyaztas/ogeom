og_features <- function(mom, alphas) {
  alphas <- as.matrix(alphas)
  P <- matrix(0, nrow(mom$Uhat), nrow(alphas))
  for (k in seq_len(nrow(alphas))) {
    a <- og_index_list(alphas[k, ])
    P[, k] <- og_wick_batch(a, mom$Uhat, mom$beta) / sqrt(prod(factorial(alphas[k, ])))
  }
  colnames(P) <- og_alpha_labels(alphas)
  P
}