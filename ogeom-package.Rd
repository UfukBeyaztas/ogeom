\name{ogeom-package}
\alias{ogeom-package}
\alias{ogeom}
\docType{package}
\title{Nonlinear Regression on Sparsely Observed Gaussian Curves}
\description{
Galerkin chaos estimators for binary regression on a Gaussian trajectory measured at a few irregular times per subject, the information operators that show which nonlinear effects of the trajectory a sampling design can identify, and the data-generating processes of the accompanying paper.
}
\details{
Functionals of the trajectory are represented through their Wiener-Ito chaos expansions. For a binary response the law of the data depends on the regression functional \eqn{g}{g} only through the conditional expectation operator \eqn{K}{K}, and under Gaussianity \eqn{K^*K = E\{\Gamma(B_D)\}}{K*K = E{Gamma(B_D)}}, where \eqn{B_D}{B_D} is a random contraction fixed by the design. The information about effects of order \eqn{q}{q} is therefore a tensor moment of \eqn{B_D}{B_D}, and components of order larger than the largest number of measurements per subject are not identified for Brownian motion and several other processes.

The main functions are \code{\link{og_fit}} and \code{\link{og_fit_constrained}} (estimation), \code{\link{og_predict}} and \code{\link{og_kernels}} (prediction and interpretation), \code{\link{og_information_profile}} (what a design reveals), \code{\link{og_process_kl}}, \code{\link{og_process_bm}} and \code{\link{og_process_from_fpca}} (law of the trajectory), \code{\link{og_bm_spectral}} and \code{\link{og_bm_augmented}} (Brownian motion measured once), and \code{\link{og_simulate}}, \code{\link{og_simulation_design}} and \code{\link{og_simulate_bm}} (simulation designs of the paper).
}
\author{Ufuk Beyaztas \email{ufuk.beyaztas@marmara.edu.tr}}
\keyword{package}
