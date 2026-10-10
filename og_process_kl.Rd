\name{og_process_kl}
\alias{og_process_kl}
\alias{og_process_bm}
\alias{og_bm_dir_kl}
\alias{og_bm_dir_const}
\alias{og_fourier_phi}
\title{Specification of the Gaussian trajectory}
\description{
Construct the process objects used by all estimation functions. A process object describes a centered Gaussian process \eqn{X(t) = W(k_t)} through its covariance function and the cross-covariances between the measurements and the latent scores \eqn{V_l = W(u_l)}{V_l = W(u_l)}, \eqn{l = 1,\dots,L}{l = 1, ..., L}.
}
\usage{
og_process_kl(mu, phi, U = diag(length(mu)))

og_process_bm(cumdir)

og_bm_dir_kl(L)

og_bm_dir_const()

og_fourier_phi(K)
}
\arguments{
  \item{mu}{positive eigenvalues \eqn{\mu_1,\dots,\mu_K}{mu_1, ..., mu_K}.}
  \item{phi}{function of a vector \code{t} returning the \code{length(t)} by \code{K} matrix of eigenfunction values.}
  \item{U}{\code{K} by \code{L} matrix with orthonormal columns defining the directions.}
  \item{cumdir}{function of \code{t} returning the \code{length(t)} by \code{L} matrix of \eqn{\int_0^t u_l(s)\,ds}{int_0^t u_l(s) ds} for orthonormal directions \eqn{u_l}{u_l} in \eqn{L^2[0,1]}{L^2[0,1]}.}
  \item{L}{number of directions.}
  \item{K}{number of eigenfunctions.}
}
\details{
\code{og_process_kl} gives a process of finite rank, \eqn{X(t) = \sum_k \mu_k^{1/2}\xi_k\phi_k(t)}{X(t) = sum_k sqrt(mu_k) xi_k phi_k(t)}, with directions given by the orthonormal columns of \code{U}; the default \code{U} makes the latent scores the standardized principal component scores. \code{og_process_bm} gives standard Brownian motion, for which the Gaussian Hilbert space is \eqn{L^2[0,1]}{L^2[0,1]}, \eqn{k_t = 1_{[0,t]}}{k_t = 1[0,t]} and \eqn{W(h) = \int h\,dX}{W(h) = int h dX}; the directions are specified through their primitives.

\code{og_bm_dir_kl}, \code{og_bm_dir_const} and \code{og_fourier_phi} supply standard choices: the Karhunen-Loeve directions of Brownian motion, the constant direction (for which \eqn{W(u) = X(1)}{W(u) = X(1)}), and Fourier eigenfunctions, the first two of which are those of the rotating process used in the paper.
}
\value{
\code{og_process_kl} and \code{og_process_bm} return a list with components \code{L} (number of directions), \code{cov} (covariance function of two time vectors) and \code{cross} (function returning the matrix of \eqn{Cov\{W(u_l), X(t_j)\}}{Cov{W(u_l), X(t_j)}}), together with the defining quantities. The other functions return functions of \code{t}.
}
\examples{
\dontrun{
## rotating process of Experiment A
phi <- function(t) cbind(sqrt(2) * cos(2 * pi * t), sqrt(2) * sin(2 * pi * t))
rot <- og_process_kl(c(0.5, 0.5), phi)
rot$cov(c(0, 0.25), c(0, 0.25))
## rank-four Fourier process of Experiment C
procC <- og_process_kl(c(1, 0.6, 0.35, 0.2), og_fourier_phi(4))
## Brownian motion with its first three Karhunen-Loeve directions
bm <- og_process_bm(og_bm_dir_kl(3))
bm$cross(c(0.2, 0.7))
}
}
