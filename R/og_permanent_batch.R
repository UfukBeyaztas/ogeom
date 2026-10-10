og_permanent_batch <- function(B, a, b, perms = NULL, method = c("auto", "ryser", "permutations")) {
  method <- match.arg(method)
  q <- length(a)
  if (length(b) != q) stop("a and b must have the same length")
  if (q == 0) return(rep(1, dim(B)[1]))
  if (method == "auto") method <- if (factorial(q) * q <= 2^q * q^2) "permutations" else "ryser"
  if (method == "permutations") {
    if (is.null(perms)) perms <- og_permutations(q)
    out <- 0
    for (k in seq_len(nrow(perms))) {
      term <- 1
      for (i in seq_len(q)) term <- term * B[, a[i], b[perms[k, i]]]
      out <- out + term
    }
    return(out)
  }
  out <- 0
  for (s in 1:(2^q - 1)) {
    cols <- which(bitwAnd(s, 2^(0:(q - 1))) > 0)
    prodv <- 1
    for (i in seq_len(q)) {
      rs <- 0
      for (j in cols) rs <- rs + B[, a[i], b[j]]
      prodv <- prodv * rs
    }
    out <- out + (-1)^length(cols) * prodv
  }
  (-1)^q * out
}
