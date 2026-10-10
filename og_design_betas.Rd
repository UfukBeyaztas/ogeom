\name{og_design_betas}
\alias{og_design_betas}
\alias{og_gram}
\alias{og_information_profile}
\title{Information operators of a sampling design}
\description{
Compute the information blocks \eqn{\mathcal I_q}{I_q} of a design, which determine which components of the regression functional the design identifies and how accurately they can be estimated.
}
\usage{
og_design_betas(process, design, sigma2 = 0, nmc = 20000, seed = NULL)

og_gram(beta, alphas, weights = NULL, method = c("auto", "ryser",
        "permutations"))

og_information_profile(process, design, sigma2 = 0, q, nmc = 20000, beta = NULL,
        tol = 1e-06, seed = NULL)
}
\arguments{
  \item{process}{process object, see \code{\link{og_process_kl}}.}
  \item{design}{design sampler: a function of \code{n} returning a list of \code{n} vectors of measurement times, as produced by \code{\link{og_design_sampler}}.}
  \item{sigma2}{measurement error variance, or a function of the times returning the error covariance matrix.}
  \item{nmc}{number of design draws.}
  \item{seed}{optional seed for the design draws; the random number state of the caller is preserved.}
  \item{beta}{array of matrices \eqn{\beta(D)}{beta(D)}, one per design (first dimension), for example the component \code{beta} of \code{\link{og_moments_all}}.}
  \item{alphas}{matrix of multi-indices of one order, see \code{\link{og_multi_indices}}.}
  \item{weights}{optional weights of the designs; equal weights by default.}
  \item{method}{algorithm for the permanents.}
  \item{q}{order of the information block.}
  \item{tol}{eigenvalues below \code{tol} are reported as weak directions.}
}
\details{
For a binary response the law of the data depends on the regression functional \eqn{g}{g} only through the conditional expectation operator \eqn{K}{K}. Under Gaussianity, \eqn{K^*K = E\{\Gamma(B_D)\}}{K*K = E{Gamma(B_D)}}, the expected second quantization of a random contraction \eqn{B_D}{B_D} fixed by the design and the measurement error. Its block of order \eqn{q}{q} is the tensor moment \eqn{E(B_D^{\otimes q})}{E(B_D^(x)q)}; in the basis \eqn{b_\alpha}{b_alpha}, its entries are expected permanents of blocks of \eqn{\beta(D)}{beta(D)}.

\code{og_design_betas} draws designs from the design law and returns their matrices \eqn{\beta(D)}{beta(D)}; \code{og_gram} averages the permanents over such matrices (or over the observed designs); \code{og_information_profile} combines the two and returns the spectrum. Zero eigenvalues correspond to components that the design cannot identify; for example, components of order \eqn{q}{q} larger than the largest number of measurements per subject are not identified for Brownian motion.
}
\value{
\code{og_design_betas} returns an \code{nmc} by \eqn{L}{L} by \eqn{L}{L} array; \code{og_gram} returns the Gram matrix \eqn{G_q}{G_q}; \code{og_information_profile} returns a list with \code{G}, its eigenvalues \code{values} and eigenvectors \code{vectors}, \code{alphas}, the indices \code{weak} of eigenvalues below \code{tol}, and \code{q}.
}
\examples{
\dontrun{
## Experiment A: equal linear information, different second-order information
for (d in c("A", "B")) {
  des <- og_simulation_design("A", design = d)
  print(round(og_information_profile(des$process, des$design, des$sigma2, q = 2)$values, 3))
}
## design A: 0.5, 0.5, 0 (the interaction V1 V2 is invisible); design B: 0.5, 0.25, 0.25
## Experiment B: second-order information grows with the share of subjects measured twice
for (p2 in c(0, 0.25, 0.5, 1)) {
  des <- og_simulation_design("B", prob2 = p2)
  ev <- og_information_profile(des$process, des$design, des$sigma2, q = 2, nmc = 5000)$values
  cat(p2, round(min(ev), 4), "\\n")
}
}
}
