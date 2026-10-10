\name{og_fb_eigen}
\alias{og_fb_eigen}
\alias{og_bessel_j0_zeros}
\alias{og_fb_phi}
\alias{og_bm_dir_fb}
\title{Fourier-Bessel eigensystem of the Brownian information operator}
\description{
Zeros of the Bessel function \eqn{J_0}{J_0}, the eigenvalues \eqn{4/j_{0,k}^2}{4 / j_{0,k}^2} and the eigenfunctions \eqn{\varphi_k(s) = J_0(j_{0,k}s^{1/2})/|J_1(j_{0,k})|}{phi_k(s) = J_0(j_{0,k} sqrt(s)) / |J_1(j_{0,k})|}.
}
\usage{
og_fb_eigen(K)

og_bessel_j0_zeros(K)

og_fb_phi(s, j)

og_bm_dir_fb(L, j = og_bessel_j0_zeros(L))
}
\arguments{
  \item{K}{number of zeros or eigenvalues.}
  \item{s}{evaluation points in \eqn{[0,1]}{[0, 1]}.}
  \item{j}{zeros of \eqn{J_0}{J_0}.}
  \item{L}{number of directions.}
}
\details{
For Brownian motion measured once, without error, at a uniformly distributed time, the first-order information operator is the integral operator with kernel \eqn{-\log\max(s,u)}{-log max(s, u)}. Its eigenvalues are \eqn{4/j_{0,k}^2}{4 / j_{0,k}^2}, which sum to one, and its eigenfunctions are the Fourier-Bessel functions above. \code{og_bm_dir_fb} gives these eigenfunctions as directions for \code{\link{og_process_bm}}.
}
\value{
\code{og_bessel_j0_zeros} returns the first \code{K} positive zeros; \code{og_fb_eigen} returns a list with \code{j} and \code{lambda}; \code{og_fb_phi} returns the matrix of eigenfunction values; \code{og_bm_dir_fb} returns a function of \code{t}.
}
\examples{
\dontrun{
e <- og_fb_eigen(5)
e$lambda
sum(og_fb_eigen(2000)$lambda)   # close to 1
s <- seq(0, 1, length.out = 401)
matplot(s, og_fb_phi(s, e$j[1:4]), type = "l", lty = 1)
}
}
