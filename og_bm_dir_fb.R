og_bm_dir_fb <- function(L, j = og_bessel_j0_zeros(L)) {
  j <- j[seq_len(L)]
  function(t) og_fb_cum(t, j)
}