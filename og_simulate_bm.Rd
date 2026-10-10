\name{og_simulate_bm}
\alias{og_simulate_bm}
\title{Brownian motion measured once at a uniform time (Experiment D)}
\description{
Generate the data of Experiment D: one noiseless measurement \eqn{W = X(T)}{W = X(T)} of standard Brownian motion at a uniform time \eqn{T}{T}, and a binary response with \eqn{P(Y = 1 \mid X) = 1/2 + \delta\sin W(u)}{P(Y = 1 | X) = 1/2 + delta sin W(u)}.
}
\usage{
og_simulate_bm(n, a = 1, b = 0, delta = 0.4)
}
\arguments{
  \item{n}{number of subjects.}
  \item{a}{coefficient of the constant part of the direction.}
  \item{b}{coefficient of \eqn{1 - s}{1 - s}.}
  \item{delta}{amplitude of the regression.}
}
\details{
The direction \eqn{u}{u} is proportional to \eqn{a + b(1-s)}{a + b(1 - s)}, and the first-order kernel of the regression functional is \eqn{\delta e^{-1/2}u}{delta exp(-1/2) u}. The kernel does not vanish at \eqn{s = 1}{s = 1} unless \eqn{a = 0}{a = 0}, which is the situation in which the augmented estimator is needed.
}
\value{
A list with \code{T}, \code{W}, \code{Y}, the latent score \code{Vu}, the values \code{g} of the regression function and the true first-order kernel \code{f1}, a function of \eqn{s}{s}.
}
\examples{
\dontrun{
set.seed(1)
d <- og_simulate_bm(2000, a = 0, b = 1)   # kernel vanishing at 1
fit <- og_bm_spectral(d$T, d$W, d$Y, og_bm_J(2000, b = 2))
s <- seq(0, 1, length.out = 201)
plot(s, d$f1(s), type = "l")
lines(s, og_bm_kernel(s, fit), col = 2)
}
}
