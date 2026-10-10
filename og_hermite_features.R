og_hermite_features <- function(V, alphas) {
  V <- as.matrix(V); alphas <- as.matrix(alphas)
  B <- matrix(1, nrow(V), nrow(alphas))
  for (k in seq_len(nrow(alphas))) {
    for (l in seq_len(ncol(alphas))) if (alphas[k, l] > 0) B[, k] <- B[, k] * og_hermite(V[, l], alphas[k, l])
    B[, k] <- B[, k] / sqrt(prod(factorial(alphas[k, ])))
  }
  colnames(B) <- og_alpha_labels(alphas)
  B
}