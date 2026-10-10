\name{og_design_sampler}
\alias{og_design_sampler}
\title{Design sampler with a random number of uniform measurement times}
\description{
Create a design sampler: the number of measurements takes the values in \code{support} with probabilities \code{probs}, and given that number the times are independent and uniform on \eqn{[0,1]}{[0, 1]}.
}
\usage{
og_design_sampler(probs, support = seq_along(probs))
}
\arguments{
  \item{probs}{probabilities of the numbers of measurements.}
  \item{support}{possible numbers of measurements.}
}
\details{
Given each number of measurements, the times have a joint law with full support, so the condition under which the null space of the information operators has its explicit diagonal description holds. With such designs, the component of order \eqn{q}{q} of the regression functional is identified for Brownian motion exactly when \eqn{q}{q} does not exceed the largest number of measurements that occurs with positive probability. The times are returned in increasing order.
}
\value{
A function of \code{n} returning a list of \code{n} increasing vectors of times.
}
\examples{
\dontrun{
set.seed(1)
draw <- og_design_sampler(c(0.4, 0.3, 0.2, 0.1))   # Experiment C
draw(3)
}
}
