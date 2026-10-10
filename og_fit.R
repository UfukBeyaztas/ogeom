og_fit <- function(data, Y, process, sigma2 = 0, Q = 2,
                   gram = c("design", "conditional", "supplied"), design = NULL,
                   nmc = 20000, beta_design = NULL, G_list = NULL,
                   select = c("ure", "full"), J = NULL, tol = 1e-10, seed = NULL) {
  gram <- match.arg(gram); select <- match.arg(select)
  n <- length(Y); L <- process$L
  if (!all(Y %in% c(0, 1))) stop("Y must be binary (0/1)")
  if (n != length(data$t) || n != length(data$w)) stop("data and Y must describe the same subjects")
  if (gram == "supplied" && length(G_list) < Q) stop("G_list must contain Gram matrices for orders 1..Q")
  mom <- og_moments_all(data, process, sigma2)
  if (gram == "design" && is.null(beta_design)) {
    if (is.null(design)) stop("gram = 'design' needs a design sampler")
    beta_design <- og_design_betas(process, design, sigma2, nmc, seed)
  }
  orders <- vector("list", Q)
  for (q in seq_len(Q)) {
    alphas <- og_multi_indices(L, q)
    Psi <- og_features(mom, alphas)
    G <- switch(gram,
                design      = og_gram(beta_design, alphas),
                conditional = og_gram(mom$beta, alphas),
                supplied    = as.matrix(G_list[[q]]))
    if (!all(dim(G) == ncol(Psi))) stop("Gram matrix of order ", q, " has the wrong dimension")
    full <- og_galerkin(Psi, Y, G, J = NULL, tol = tol)
    ure <- og_ure(full$eta, full$values[seq_len(full$npos)], n)
    Jq <- if (!is.null(J)) J[min(q, length(J))] else if (select == "full") full$npos else which.min(ure) - 1
    fq <- og_galerkin(Psi, Y, G, J = Jq, tol = tol)
    fq$alphas <- alphas; fq$G <- G; fq$ure <- ure
    fq$covariance_scope <- if (gram == "conditional" || (is.null(J) && select == "ure"))
      "Diagnostic only: A depends on the fitted sample; not a valid covariance for inference" else
        "Fixed A independent of the observed sample: sandwich covariance estimates sampling variation, not bias"
    orders[[q]] <- fq
  }
  criterion <- if (!is.null(J)) "fixed truncation" else if (select == "full") "all positive eigen-directions" else
    switch(gram,
           supplied    = "risk criterion with the supplied eigenpairs (exact for each fixed J if they are the population ones, prop:ure)",
           design      = "approximate spectral risk criterion (Monte Carlo eigenpairs)",
           conditional = "heuristic tuning criterion (design-conditional Gram)")
  structure(list(orders = orders, Ybar = mean(Y), Q = Q, process = process,
                 sigma2 = sigma2, gram = gram, criterion = criterion, n = n), class = "og_fit")
}