og_wick_batch <- function(a, Uhat, beta) {
  p <- length(a)
  if (p == 0) return(rep(1, nrow(Uhat)))
  rest <- a[-p]; last <- a[p]
  val <- Uhat[, last] * og_wick_batch(rest, Uhat, beta)
  for (i in seq_len(p - 1)) val <- val - beta[, rest[i], last] * og_wick_batch(rest[-i], Uhat, beta)
  val
}
