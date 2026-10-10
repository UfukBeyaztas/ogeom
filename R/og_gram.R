og_gram <- function(beta, alphas, weights = NULL, method = c("auto", "ryser", "permutations")) {
  method <- match.arg(method)
  alphas <- as.matrix(alphas)
  N <- dim(beta)[1]; p <- nrow(alphas); q <- sum(alphas[1, ])
  if (any(rowSums(alphas) != q)) stop("all multi-indices must have the same order")
  if (is.null(weights)) weights <- rep(1 / N, N)
  if (length(weights) != N || any(weights < 0) || abs(sum(weights) - 1) > 1e-8)
    stop("weights must be nonnegative and sum to one")
  perms <- if (method != "ryser" && q <= 6) og_permutations(q) else NULL
  idx <- lapply(seq_len(p), function(k) og_index_list(alphas[k, ]))
  fac <- sqrt(apply(alphas, 1, function(a) prod(factorial(a))))
  G <- matrix(0, p, p)
  for (k in seq_len(p)) for (k2 in k:p) {
    pv <- og_permanent_batch(beta, idx[[k]], idx[[k2]], perms, method)
    G[k, k2] <- G[k2, k] <- sum(weights * pv) / (fac[k] * fac[k2])
  }
  dimnames(G) <- list(og_alpha_labels(alphas), og_alpha_labels(alphas))
  G
}
