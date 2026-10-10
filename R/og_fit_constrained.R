og_fit_constrained <- function(data, Y, process, sigma2 = 0, Q = 2,
                               gram = c("design", "conditional", "supplied"), design = NULL,
                               nmc = 20000, beta_design = NULL, G_list = NULL,
                               bound = 0.25, tol = 1e-10, seed = NULL) {
  gram <- match.arg(gram)
  if (gram == "design" && is.null(beta_design)) {
    if (is.null(design)) stop("gram = 'design' needs a design sampler")
    beta_design <- og_design_betas(process, design, sigma2, nmc, seed)
  }
  full <- og_fit(data, Y, process, sigma2, Q, gram = gram, beta_design = beta_design,
                 G_list = G_list, select = "full", tol = tol)
  J <- vapply(full$orders, function(o) {
    if (o$npos == 0) return(0)
    ure <- og_ure(o$eta, o$values[seq_len(o$npos)], length(Y))
    ok  <- c(TRUE, cumsum(o$eta^2) <= bound)
    as.numeric(which.min(ifelse(ok, ure, Inf)) - 1)
  }, numeric(1))
  fit <- og_fit(data, Y, process, sigma2, Q, gram = gram, beta_design = beta_design,
                G_list = G_list, J = J, tol = tol)
  fit$criterion <- sprintf("risk criterion of prop:ure restricted to ||g_q||^2 <= %g", bound)
  fit
}
