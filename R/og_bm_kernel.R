og_bm_kernel <- function(s, fit) {
  f <- if (length(fit$theta) > 0) drop(og_fb_phi(s, fit$j) %*% fit$theta) else rep(0, length(s))
  if (!is.null(fit$theta_e)) {
    proj <- if (length(fit$j) > 0) drop(og_fb_phi(s, fit$j) %*% fit$one_k) else 0
    f <- f + fit$theta_e * (1 - proj) / fit$rJ
  }
  f
}
