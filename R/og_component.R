og_component <- function(fit, q, V) drop(og_hermite_features(V, fit$orders[[q]]$alphas) %*% fit$orders[[q]]$theta)
