print.og_fit <- function(x, ...) {
  cat("Galerkin chaos estimator:", x$n, "subjects, orders 1 to", x$Q, "| Gram:", x$gram, "\n")
  cat("  truncation:", x$criterion, "\n")
  for (q in seq_len(x$Q)) {
    o <- x$orders[[q]]
    cat(sprintf("  order %d: %d coefficients, %d positive eigenvalues, J = %d retained\n",
                q, length(o$theta), o$npos, o$J))
  }
  invisible(x)
}