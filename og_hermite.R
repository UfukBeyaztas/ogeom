og_hermite <- function(x, m) {
  if (m == 0) return(rep(1, length(x)))
  h0 <- rep(1, length(x)); h1 <- x
  if (m == 1) return(h1)
  for (k in 1:(m - 1)) {
    h2 <- x * h1 - k * h0
    h0 <- h1; h1 <- h2
  }
  h1
}