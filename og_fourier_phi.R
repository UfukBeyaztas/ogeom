og_fourier_phi <- function(K) {
  force(K)
  function(t) {
    out <- matrix(0, length(t), K)
    for (k in seq_len(K)) {
      fr <- 2 * pi * ceiling(k / 2)
      out[, k] <- if (k %% 2 == 1) sqrt(2) * sin(fr * t) else sqrt(2) * cos(fr * t)
    }
    out
  }
}