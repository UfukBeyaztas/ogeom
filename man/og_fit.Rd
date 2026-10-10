\name{og_fit}
\alias{og_fit}
\alias{print.og_fit}
\title{Galerkin chaos estimator}
\description{
Estimate the chaos components of a binary regression functional of a sparsely observed Gaussian trajectory.
}
\usage{
og_fit(data, Y, process, sigma2 = 0, Q = 2, gram = c("design",
        "conditional", "supplied"), design = NULL, nmc = 20000, beta_design = NULL,
        G_list = NULL, select = c("ure", "full"), J = NULL, tol = 1e-10,
        seed = NULL)

\method{print}{og_fit}(x, \dots)
}
\arguments{
  \item{data}{list with components \code{t} and \code{w}: lists of measurement times and centered measurements, one element per subject.}
  \item{Y}{binary responses (0 or 1).}
  \item{process}{process object, see \code{\link{og_process_kl}}.}
  \item{sigma2}{measurement error variance, or a function of the times returning the error covariance matrix.}
  \item{Q}{highest chaos order to estimate.}
  \item{gram}{source of the Gram matrices: \code{"design"} uses the known design law through design draws independent of the data; \code{"conditional"} averages the exact conditional Gram matrices of the observed designs; \code{"supplied"} uses \code{G_list}.}
  \item{design}{design sampler, needed for \code{gram = "design"} unless \code{beta_design} is given.}
  \item{nmc}{number of design draws for \code{gram = "design"}.}
  \item{beta_design}{optional array of design matrices \eqn{\beta(D)}{beta(D)} from \code{\link{og_design_betas}}, reused across fits.}
  \item{G_list}{list of Gram matrices of orders \code{1:Q}, for \code{gram = "supplied"}.}
  \item{select}{truncation rule: \code{"ure"} minimizes the risk criterion; \code{"full"} uses all positive eigen-directions.}
  \item{J}{optional vector of fixed truncation levels, overriding \code{select}.}
  \item{tol}{relative tolerance for positive eigenvalues.}
  \item{seed}{optional seed for the design draws.}
  \item{x}{object of class \code{og_fit}.}
  \item{\dots}{unused.}
}
\details{
Write \eqn{g = \sum_q g_q}{g = sum_q g_q} for the Wiener-Ito chaos decomposition of the regression functional and expand \eqn{g_q = \sum_{|\alpha| = q} \theta_\alpha b_\alpha(V)}{g_q = sum over |alpha| = q of theta_alpha b_alpha(V)}. For a binary response, \eqn{E\{Y\psi_\alpha(O)\} = \langle Kg, Kb_\alpha\rangle}{E{Y psi_alpha(O)} = <Kg, K b_alpha>}, and since \eqn{K^*K}{K*K} preserves the chaos decomposition, the coefficients of order \eqn{q}{q} solve \eqn{G_q\theta_q = c_q}{G_q theta_q = c_q}, free of the components of other orders.

The estimator replaces \eqn{c_q}{c_q} by its empirical version and inverts \eqn{G_q}{G_q} on its leading \eqn{J}{J} eigen-directions. Its mean squared error has an exact finite-sample expression, and the truncation level is chosen for each order by an unbiased estimate of the risk up to a constant. Prediction uses \eqn{P(Y = 1 \mid O) = \bar Y + \sum_q \psi_q(O)^\top\hat\theta_q}{P(Y = 1 | O) = mean(Y) + sum_q psi_q(O)' theta_q}. Components outside the range of the information operators, for instance those above the identification threshold, are not estimable by any method.
}
\value{
An object of class \code{og_fit}: a list with \code{orders} (one list per order with the coefficients \code{theta}, the truncation level \code{J}, the eigenvalues \code{values} and eigenvectors \code{vectors} of the Gram matrix \code{G}, the multi-indices \code{alphas} and the risk criterion \code{ure}), the mean response \code{Ybar}, \code{Q}, \code{process}, \code{sigma2}, \code{gram}, the truncation \code{criterion} and \code{n}.
}
\examples{
\dontrun{
## Experiment C, interaction scenario
set.seed(1)
des   <- og_simulation_design("C", scenario = "interaction")
train <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
test  <- og_simulate(2000, des$process, des$design, des$sigma2, des$g, des$r)
fit <- og_fit(train, train$Y, des$process, des$sigma2, Q = 2, gram = "design", design = des$design)
fit
p <- og_predict(fit, test)
sqrt(mean((p - test$Kg)^2))   # distance to the Bayes predictor P(Y = 1 | O)
}
}
