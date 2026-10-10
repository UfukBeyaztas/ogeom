og_design_sampler <- function(probs, support = seq_along(probs)) {
  if (any(probs < 0) || abs(sum(probs) - 1) > 1e-8) stop("probs must be nonnegative and sum to one")
  force(support)
  function(n) {
    M <- support[sample.int(length(support), n, replace = TRUE, prob = probs)]
    lapply(M, function(m) sort(runif(m)))
  }
}