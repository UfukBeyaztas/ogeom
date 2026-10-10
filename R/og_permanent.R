og_permanent <- function(M) {
  M <- as.matrix(M); q <- nrow(M)
  if (q == 0) return(1)
  total <- 0
  for (s in 1:(2^q - 1)) {
    cols <- which(bitwAnd(s, 2^(0:(q - 1))) > 0)
    total <- total + (-1)^length(cols) * prod(rowSums(M[, cols, drop = FALSE]))
  }
  (-1)^q * total
}
