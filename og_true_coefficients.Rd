\name{og_true_coefficients}
\alias{og_true_coefficients}
\alias{og_oracle_Kg}
\title{True chaos coefficients and the Bayes predictor}
\description{
Compute, for simulated data, the true chaos coefficients \eqn{\theta_\alpha = E\{g(V)b_\alpha(V)\}}{theta_alpha = E{g(V) b_alpha(V)}} and the Bayes predictor \eqn{P(Y = 1 \mid O) = E\{g(V) \mid O\}}{P(Y = 1 | O) = E{g(V) | O}}.
}
\usage{
og_true_coefficients(gfun, r, alphas, nq = 30)

og_oracle_Kg(mom, gfun, r, nq = 20)
}
\arguments{
  \item{gfun}{regression function of a matrix of latent scores.}
  \item{r}{number of latent scores on which \code{gfun} depends.}
  \item{alphas}{matrix of multi-indices, see \code{\link{og_multi_indices}}.}
  \item{nq}{number of quadrature nodes per coordinate.}
  \item{mom}{list returned by \code{\link{og_moments_all}}.}
}
\details{
Both quantities are computed by tensor Gauss-Hermite quadrature. Given \eqn{O}{O}, the latent scores are normal with mean \eqn{\hat U}{Uhat} and covariance \eqn{I - \beta(D)}{I - beta(D)}, so the Bayes predictor is a Gaussian integral of \eqn{g}{g}. They serve to evaluate estimators in simulations.
}
\value{
\code{og_true_coefficients} returns the named vector of coefficients; \code{og_oracle_Kg} returns the vector of probabilities.
}
\examples{
\dontrun{
sc <- og_scenario("interaction")
og_true_coefficients(sc$g, 2, og_multi_indices(2, 2))   # V1*V2: 0.45 exp(-1)
}
}
