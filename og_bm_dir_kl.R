og_bm_dir_kl <- function(L) {
  force(L)
  function(t) {
    w <- (seq_len(L) - 0.5) * pi
    sqrt(2) * sin(outer(t, w)) / rep(w, each = length(t))
  }
}