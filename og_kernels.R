og_kernels <- function(fit, grid = seq(0, 1, length.out = 101)) {
  pr <- fit$process
  if (!identical(pr$type, "kl")) stop("og_kernels needs a process created by og_process_kl")
  S   <- matrix(pr$phi(grid), nrow = length(grid), ncol = pr$K) %*% diag(1 / sqrt(pr$mu), pr$K)
  out <- list(grid = grid, beta1 = drop(S %*% (pr$U %*% fit$orders[[1]]$theta)))
  if (fit$Q >= 2) {
    al <- fit$orders[[2]]$alphas; th <- fit$orders[[2]]$theta; B <- matrix(0, pr$L, pr$L)
    for (r in seq_len(nrow(al))) {
      ix <- which(al[r, ] > 0)
      if (length(ix) == 1) B[ix, ix] <- th[r] / sqrt(2) else B[ix[1], ix[2]] <- B[ix[2], ix[1]] <- th[r] / 2
    }
    out$beta2 <- S %*% pr$U %*% B %*% t(pr$U) %*% t(S)
    out$B <- B
  }
  out
}