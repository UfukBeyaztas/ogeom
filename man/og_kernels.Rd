\name{og_kernels}
\alias{og_kernels}
\title{First- and second-order kernels of a fitted regression functional}
\description{
Express the fitted first- and second-order components through kernels on the time interval.
}
\usage{
og_kernels(fit, grid = seq(0, 1, length.out = 101))
}
\arguments{
  \item{fit}{object of class \code{og_fit} fitted with a process from \code{\link{og_process_kl}}.}
  \item{grid}{evaluation grid in \eqn{[0,1]}{[0, 1]}.}
}
\details{
For a process of finite rank, the fitted components satisfy \eqn{g_1 = \int\beta_1(t)X(t)\,dt}{g_1 = int beta_1(t) X(t) dt} and \eqn{g_2 = \int\int\beta_2(s,t)\{X(s)X(t) - C(s,t)\}\,ds\,dt}{g_2 = int int beta_2(s, t) {X(s) X(t) - C(s, t)} ds dt}. The kernels show which parts of the trajectory, and which pairs of times, drive the risk. The eigenfunctions are assumed orthonormal in \eqn{L^2[0,1]}{L^2[0, 1]}.
}
\value{
A list with \code{grid}, the vector \code{beta1}, and, when \code{fit$Q >= 2}, the matrix \code{beta2} and the symmetric coefficient matrix \code{B}.
}
\examples{
\dontrun{
## Experiment C, interaction scenario
set.seed(1)
des   <- og_simulation_design("C", scenario = "interaction")
train <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
fit <- og_fit(train, train$Y, des$process, des$sigma2, Q = 2, design = des$design)
k <- og_kernels(fit)
image(k$grid, k$grid, k$beta2, xlab = "s", ylab = "t")
}
}
