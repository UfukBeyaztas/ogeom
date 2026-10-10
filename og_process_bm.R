og_process_bm <- function(cumdir) {
  L <- ncol(as.matrix(cumdir(0.5)))
  list(type = "bm", L = L, cumdir = cumdir,
       cov   = function(s, t) outer(s, t, pmin),
       cross = function(t) matrix(cumdir(t), nrow = length(t), ncol = L))
}