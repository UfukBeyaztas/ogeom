\name{og_bm_spectral}
\alias{og_bm_spectral}
\alias{og_bm_augmented}
\alias{og_bm_kernel}
\alias{og_bm_J}
\title{Fourier-Bessel estimators for Brownian motion measured once}
\description{
Estimate the first-order kernel of a binary regression functional of Brownian motion observed once, without error, at a uniformly distributed time.
}
\usage{
og_bm_spectral(T, W, Y, J, j = og_bessel_j0_zeros(J))

og_bm_augmented(T, W, Y, J, j_all = og_bessel_j0_zeros(J + 2000))

og_bm_kernel(s, fit)

og_bm_J(n, b)
}
\arguments{
  \item{T}{measurement times.}
  \item{W}{measurements \eqn{X(T)}{X(T)}.}
  \item{Y}{binary responses.}
  \item{J}{number of Fourier-Bessel terms.}
  \item{j}{the first \code{J} positive zeros of \eqn{J_0}{J_0}.}
  \item{j_all}{positive zeros of \eqn{J_0}{J_0}, at least \code{J} of them; the remaining ones are used for the tails of the augmented basis.}
  \item{s}{evaluation points in \eqn{[0,1]}{[0, 1]}.}
  \item{fit}{list returned by \code{og_bm_spectral} or \code{og_bm_augmented}.}
  \item{n}{sample size.}
  \item{b}{smoothness index of the ellipsoid.}
}
\details{
For this design the first-order information operator has kernel \eqn{-\log\max(s,u)}{-log max(s, u)}, eigenvalues \eqn{\lambda_k = 4/j_{0,k}^2}{lambda_k = 4 / j_{0,k}^2} and Fourier-Bessel eigenfunctions, see \code{\link{og_fb_eigen}}. The spectral estimator truncates the expansion at \eqn{J}{J} terms and attains the minimax rate \eqn{n^{-2b/(2b+3)}}{n^(-2b/(2b+3))} over Fourier-Bessel ellipsoids with the truncation level given by \code{og_bm_J}.

For smoothness \eqn{b \ge 1/2}{b >= 1/2} the ellipsoids impose the boundary condition \eqn{f_1(1) = 0}{f_1(1) = 0}. When the kernel does not vanish at 1 the spectral estimator cannot do better than \eqn{n^{-1/4}}{n^(-1/4)}; the augmented estimator adds the constant function to the sieve and restores the rate for twice differentiable kernels.
}
\value{
\code{og_bm_spectral} and \code{og_bm_augmented} return lists with the coefficients \code{theta} (and \code{theta_e}, the coefficient of the added function, for the augmented estimator), the zeros \code{j} and eigenvalues \code{lambda}. \code{og_bm_kernel} returns the estimated kernel at \code{s}; \code{og_bm_J} returns \eqn{\lceil n^{1/(2b+3)}\rceil}{ceiling(n^(1/(2b+3)))}.
}
\examples{
\dontrun{
## Experiment D, constant kernel (which does not vanish at 1)
set.seed(1)
d <- og_simulate_bm(2000, a = 1, b = 0)
J <- og_bm_J(2000, b = 2)
fit_s <- og_bm_spectral(d$T, d$W, d$Y, J)
fit_a <- og_bm_augmented(d$T, d$W, d$Y, J)
s <- seq(0, 1, length.out = 201)
plot(s, d$f1(s), type = "l", ylim = c(-0.1, 0.4), ylab = "first-order kernel")
lines(s, og_bm_kernel(s, fit_s), col = 2)
lines(s, og_bm_kernel(s, fit_a), col = 4)
}
}
