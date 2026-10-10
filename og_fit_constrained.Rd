\name{og_fit_constrained}
\alias{og_fit_constrained}
\title{Galerkin chaos estimator with bound-constrained truncation}
\description{
Fit the Galerkin chaos estimator with truncation levels restricted to fitted components that satisfy \eqn{\|\hat g_q\|^2 \le}{||g_q||^2 <=} \code{bound}.
}
\usage{
og_fit_constrained(data, Y, process, sigma2 = 0, Q = 2, gram = c("design",
        "conditional", "supplied"), design = NULL, nmc = 20000, beta_design = NULL,
        G_list = NULL, bound = 0.25, tol = 1e-10, seed = NULL)
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
  \item{bound}{upper bound for \eqn{\|\hat g_q\|^2}{||g_q||^2}.}
  \item{tol}{relative tolerance for positive eigenvalues.}
  \item{seed}{optional seed for the design draws.}
}
\details{
Every regression function with values in \eqn{[0,1]}{[0, 1]} satisfies \eqn{\|g_q\|^2 \le Var\, g \le 1/4}{||g_q||^2 <= Var g <= 1/4}, so the default bound adds no assumption. For each order, the truncation level minimizes the risk criterion among the levels whose fitted component satisfies the bound. The restriction protects against eigen-directions dominated by misspecification when the law of the trajectory has been estimated, as in real-data applications.
}
\value{
An object of class \code{og_fit}, see \code{\link{og_fit}}.
}
\examples{
\dontrun{
## Experiment C, interaction scenario
set.seed(1)
des   <- og_simulation_design("C", scenario = "interaction")
train <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
test  <- og_simulate(2000, des$process, des$design, des$sigma2, des$g, des$r)
fit <- og_fit_constrained(train, train$Y, des$process, des$sigma2, Q = 2, design = des$design)
sqrt(mean((og_predict(fit, test) - test$Kg)^2))
}
}
