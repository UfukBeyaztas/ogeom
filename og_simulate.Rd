\name{og_simulate}
\alias{og_simulate}
\title{Simulate sparse functional data with a binary response}
\description{
Generate measurement times, noisy measurements of a Gaussian trajectory and binary responses from the exact joint Gaussian law, without discretizing the trajectory.
}
\usage{
og_simulate(n, process, design, sigma2 = 0, gfun, r, oracle = TRUE,
        nq = 20)
}
\arguments{
  \item{n}{number of subjects.}
  \item{process}{process object, see \code{\link{og_process_kl}}.}
  \item{design}{design sampler, see \code{\link{og_design_sampler}}.}
  \item{sigma2}{measurement error variance, or a function of the times returning the error covariance matrix.}
  \item{gfun}{regression function of an \eqn{n \times r}{n x r} matrix of latent scores, with values in \eqn{[0,1]}{[0, 1]}.}
  \item{r}{number of latent scores on which \code{gfun} depends.}
  \item{oracle}{if \code{TRUE}, also compute \eqn{P(Y = 1 \mid O)}{P(Y = 1 | O)}.}
  \item{nq}{number of quadrature nodes per coordinate.}
}
\details{
The latent scores \eqn{V = (W(u_1),\dots,W(u_L))}{V = (W(u_1), ..., W(u_L))} are drawn from \eqn{N(0, I_L)}{N(0, I_L)}, and given \eqn{V}{V} the values of the trajectory at the design times are drawn from their conditional Gaussian law. Measurement errors are added, and \eqn{Y}{Y} is drawn with \eqn{P(Y = 1 \mid X) = g(V_1,\dots,V_r)}{P(Y = 1 | X) = g(V_1, ..., V_r)}. Optionally, the Bayes predictor \eqn{P(Y = 1 \mid O)}{P(Y = 1 | O)} is computed by Gauss-Hermite quadrature for evaluation. The settings of the paper's simulations are available through \code{\link{og_simulation_design}}.
}
\value{
A list with \code{t} and \code{w} (lists of times and measurements), \code{Y}, the latent scores \code{V}, \code{g} (the values of the regression function) and, if \code{oracle = TRUE}, \code{Kg}.
}
\examples{
\dontrun{
set.seed(1)
des <- og_simulation_design("C", scenario = "interaction")
dat <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
table(lengths(dat$t))   # numbers of measurements per subject
mean(dat$Y)
}
}
