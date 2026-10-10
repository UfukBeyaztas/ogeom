\name{og_process_from_fpca}
\alias{og_process_from_fpca}
\alias{og_center_data}
\title{Estimated law of the trajectory for real data}
\description{
Turn a principal component analysis through conditional expectation (PACE), as computed by \code{FPCA} in package \pkg{fdapace}, into a process object, and center the measurements with the estimated mean.
}
\usage{
og_process_from_fpca(fp, K = fp$selectK, sigma2_min = 0)

og_center_data(data, mu_fun)
}
\arguments{
  \item{fp}{object returned by \code{fdapace::FPCA}.}
  \item{K}{number of components; by default the number selected by \code{FPCA}.}
  \item{sigma2_min}{lower bound for the error variance.}
  \item{data}{list with components \code{t} and \code{w}, the lists of measurement times and measurements of the subjects.}
  \item{mu_fun}{mean function of \code{t}, for example the component \code{mu_fun} returned by \code{og_process_from_fpca}.}
}
\details{
The theory treats the law of \eqn{X} as known. With real data, the mean, eigenvalues, eigenfunctions and error variance are estimated from the pooled measurements, without the responses, and then treated as known. The signs of the eigenfunctions are fixed so that each has a positive mean. A positive floor for the error variance is advisable when the PACE estimate is close to zero, because a zero error variance makes the conditional features of nearly coincident measurement times unstable.
}
\value{
\code{og_process_from_fpca} returns a list with components \code{process} (a process object), \code{mu_fun}, \code{sigma2}, \code{lambda}, \code{phi_fun} and \code{sign}. \code{og_center_data} returns \code{data} with centered measurements.
}
\examples{
\dontrun{
library(fdapace)
set.seed(1)
des <- og_simulation_design("C", scenario = "usual")
dat <- og_simulate(500, des$process, des$design, des$sigma2, des$g, des$r, oracle = FALSE)
fp  <- FPCA(dat$w, dat$t, list(dataType = "Sparse", error = TRUE, FVEthreshold = 0.99))
law <- og_process_from_fpca(fp, sigma2_min = 0.01)
cdat <- og_center_data(dat, law$mu_fun)
fit <- og_fit_constrained(cdat, dat$Y, law$process, law$sigma2, Q = 2, gram = "conditional")
fit
}
}
