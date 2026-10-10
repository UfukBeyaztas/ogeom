\name{og_multi_indices}
\alias{og_multi_indices}
\alias{og_hermite_features}
\title{Multi-indices and Hermite features}
\description{
List the multi-indices of a given order and evaluate the normalized Hermite features \eqn{b_\alpha}{b_alpha}.
}
\usage{
og_multi_indices(L, q)

og_hermite_features(V, alphas)
}
\arguments{
  \item{L}{number of latent scores.}
  \item{q}{order of the chaos.}
  \item{V}{matrix with one row per point and \code{ncol(alphas)} columns.}
  \item{alphas}{matrix of multi-indices, one per row.}
}
\details{
The multi-indices \eqn{\alpha}{alpha} with \eqn{|\alpha| = q}{|alpha| = q} index an orthonormal basis \eqn{b_\alpha(V) = \prod_l He_{\alpha_l}(V_l)/(\alpha!)^{1/2}}{b_alpha(V) = prod_l He_alpha_l(V_l) / sqrt(alpha!)} of the chaos of order \eqn{q}{q} spanned by \eqn{L}{L} latent scores, where \eqn{He}{He} are the probabilists' Hermite polynomials. There are \code{choose(L + q - 1, q)} of them.
}
\value{
\code{og_multi_indices} returns a matrix with one multi-index per row. \code{og_hermite_features} returns the matrix of features \eqn{b_\alpha(V)}{b_alpha(V)}.
}
\examples{
\dontrun{
og_multi_indices(2, 2)
og_hermite_features(matrix(c(0.7, -0.3), 1), og_multi_indices(2, 2))
}
}
