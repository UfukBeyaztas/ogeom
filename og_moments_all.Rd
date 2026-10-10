\name{og_moments_all}
\alias{og_moments_all}
\alias{og_features}
\title{Conditional moments and exact conditional features}
\description{
Compute, for each subject, the conditional means \eqn{\hat U = E(V \mid O)}{Uhat = E(V | O)} and the matrix \eqn{\beta(D)}{beta(D)}, and from them the exact conditional features \eqn{\psi_\alpha(O) = E\{b_\alpha(V) \mid O\}}{psi_alpha(O) = E{b_alpha(V) | O}}.
}
\usage{
og_moments_all(data, process, sigma2 = 0, tol = 1e-10)

og_features(mom, alphas)
}
\arguments{
  \item{data}{list with components \code{t} and \code{w}: lists of measurement times and centered measurements, one element per subject.}
  \item{process}{process object, see \code{\link{og_process_kl}}.}
  \item{sigma2}{measurement error variance, or a function of the times returning the error covariance matrix.}
  \item{tol}{relative tolerance for the generalized inverse.}
  \item{mom}{list returned by \code{og_moments_all}.}
  \item{alphas}{matrix of multi-indices, one per row, as returned by \code{\link{og_multi_indices}}.}
}
\details{
For a design with times \eqn{t_1,\dots,t_r}{t_1, ..., t_r} and measurements \eqn{W_j = X(t_j) + \eta_j}{W_j = X(t_j) + eta_j}, the latent scores are, given \eqn{O}, normal with mean \eqn{\hat U}{Uhat} and covariance \eqn{I - \beta(D)}{I - beta(D)}, where \eqn{\beta_{jl}(D) = \langle u_j, B_D u_l\rangle}{beta_jl(D) = <u_j, B_D u_l>}. The features \eqn{b_\alpha(V) = \prod_l He_{\alpha_l}(V_l)/(\alpha!)^{1/2}}{b_alpha(V) = prod_l He_alpha_l(V_l) / sqrt(alpha!)} form an orthonormal basis of the chaos of order \eqn{|\alpha|}{|alpha|}, and their conditional expectations are Wick polynomials of \eqn{\hat U}{Uhat} with covariance \eqn{\beta(D)}{beta(D)}. These features are exact: no plug-in approximation \eqn{b_\alpha(\hat U)}{b_alpha(Uhat)} is involved.
}
\value{
\code{og_moments_all} returns a list with \code{Uhat} (\eqn{n \times L}{n x L}) and \code{beta} (\eqn{n \times L \times L}{n x L x L} array). \code{og_features} returns the \eqn{n}{n} by \code{nrow(alphas)} matrix of features.
}
\examples{
\dontrun{
set.seed(1)
des <- og_simulation_design("C")
dat <- og_simulate(200, des$process, des$design, des$sigma2, des$g, des$r, oracle = FALSE)
mom <- og_moments_all(dat, des$process, des$sigma2)
psi <- og_features(mom, og_multi_indices(4, 2))
dim(psi)
}
}
