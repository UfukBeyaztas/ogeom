og_information_profile <- function(process, design, sigma2 = 0, q, nmc = 20000,
                                   beta = NULL, tol = 1e-6, seed = NULL) {
  if (is.null(beta)) beta <- og_design_betas(process, design, sigma2, nmc, seed)
  alphas <- og_multi_indices(process$L, q)
  G <- og_gram(beta, alphas)
  e <- eigen(G, symmetric = TRUE)
  list(G = G, values = e$values, vectors = e$vectors, alphas = alphas,
       weak = which(e$values < tol), q = q)
}
