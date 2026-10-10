og_alpha_labels <- function(alphas) {
  alphas <- as.matrix(alphas)
  apply(alphas, 1, function(a) {
    s <- character(0)
    for (l in seq_along(a)) {
      if (a[l] == 1) s <- c(s, paste0("V", l))
      if (a[l] > 1) s <- c(s, paste0("H", a[l], "(V", l, ")"))
    }
    if (length(s) == 0) "1" else paste(s, collapse = "*")
  })
}
