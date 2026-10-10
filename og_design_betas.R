og_design_betas <- function(process, design, sigma2 = 0, nmc = 20000, seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
    set.seed(seed)
  }
  times <- design(nmc)
  L <- process$L
  beta <- array(0, c(nmc, L, L))
  for (s in seq_len(nmc)) {
    tt <- times[[s]]
    if (length(tt) > 0) beta[s, , ] <- og_conditional_moments(tt, numeric(length(tt)), process, sigma2)$beta
  }
  beta
}