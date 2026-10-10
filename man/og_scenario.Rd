\name{og_scenario}
\alias{og_scenario}
\title{Regression functions of the simulations}
\description{
Regression functions used in the paper's simulations, as functions of the latent scores.
}
\usage{
og_scenario(scenario = c("interaction", "usual"), delta = 0.45,
        slope = c(1.5, -1, 0.5))
}
\arguments{
  \item{scenario}{\code{"interaction"} or \code{"usual"}.}
  \item{delta}{amplitude of the interaction.}
  \item{slope}{slopes of the logistic index.}
}
\details{
\code{"interaction"}: \eqn{g = 1/2 + \delta\sin(V_1)\sin(V_2)}{g = 1/2 + delta sin(V_1) sin(V_2)}. It has no linear component, and its component of order two is \eqn{\delta e^{-1}V_1V_2}{delta exp(-1) V_1 V_2}, so it is invisible to linear methods and identified only when subjects with at least two measurements occur. \code{"usual"}: the single-index logistic model \eqn{g = \mathrm{logit}^{-1}(\sum_l s_l V_l)}{g = plogis(sum_l s_l V_l)}, dominated by its linear component.
}
\value{
A list with \code{g}, a function of a matrix of latent scores, and \code{r}, the number of scores it uses.
}
\examples{
\dontrun{
sc <- og_scenario("interaction")
sc$g(matrix(c(1, 1), 1))
}
}
