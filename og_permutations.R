og_permutations <- function(q) {
  if (q == 1) return(matrix(1, 1, 1))
  P <- og_permutations(q - 1)
  out <- NULL
  for (i in seq_len(nrow(P))) for (pos in 0:(q - 1)) out <- rbind(out, append(P[i, ], q, after = pos))
  out
}