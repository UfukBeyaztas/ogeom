\name{og_predict}
\alias{og_predict}
\alias{og_component}
\title{Prediction and fitted chaos components}
\description{
Predict \eqn{P(Y = 1 \mid O)}{P(Y = 1 | O)} for new subjects, and evaluate a fitted chaos component as a function of the latent scores.
}
\usage{
og_predict(fit, newdata, clip = TRUE)

og_component(fit, q, V)
}
\arguments{
  \item{fit}{object of class \code{og_fit}.}
  \item{newdata}{list with components \code{t} and \code{w} for the new subjects.}
  \item{clip}{if \code{TRUE}, predictions are clipped to \eqn{[0,1]}{[0, 1]}.}
  \item{q}{order of the component.}
  \item{V}{matrix of latent scores, one row per point.}
}
\details{
The prediction is \eqn{\bar Y + \sum_q \psi_q(O)^\top\hat\theta_q}{mean(Y) + sum_q psi_q(O)' theta_q}, where \eqn{\psi_q(O)}{psi_q(O)} are the exact conditional features of the new subjects. \code{og_component} returns \eqn{\hat g_q(V) = \sum_\alpha \hat\theta_\alpha b_\alpha(V)}{g_q(V) = sum_alpha theta_alpha b_alpha(V)}.
}
\value{
A vector of predicted probabilities, or of component values.
}
\examples{
\dontrun{
## Experiment C, interaction scenario
set.seed(1)
des   <- og_simulation_design("C", scenario = "interaction")
train <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
test  <- og_simulate(2000, des$process, des$design, des$sigma2, des$g, des$r)
fit <- og_fit(train, train$Y, des$process, des$sigma2, Q = 2, design = des$design)
head(og_predict(fit, test))
og_component(fit, 2, cbind(c(1, -1), c(1, 1), 0, 0))   # close to 0.45 exp(-1) V1 V2
}
}
