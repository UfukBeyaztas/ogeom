og_scenario <- function(scenario = c("interaction", "usual"), delta = 0.45,
                        slope = c(1.5, -1, 0.5)) {
  scenario <- match.arg(scenario)
  if (scenario == "interaction") {
    list(g = function(V) 0.5 + delta * sin(V[, 1]) * sin(V[, 2]), r = 2)
  } else {
    list(g = function(V) plogis(drop(as.matrix(V)[, seq_along(slope), drop = FALSE] %*% slope)),
         r = length(slope))
  }
}
