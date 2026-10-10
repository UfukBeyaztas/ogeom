og_multi_indices <- function(L, q) {
  if (L == 1) return(matrix(q, 1, 1))
  out <- NULL
  for (a1 in q:0) out <- rbind(out, cbind(a1, og_multi_indices(L - 1, q - a1)))
  unname(out)
}
