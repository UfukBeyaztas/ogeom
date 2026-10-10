# Brownian motion measured once at a uniform time (Experiment D)
#
# Each subject has one noiseless measurement W = X(T) with T uniform on [0, 1],
# and P(Y = 1 | X) = 1/2 + delta sin W(u), where u is proportional to
# a + b(1 - s). The first-order kernel is delta exp(-1/2) u.
og_simulate_bm <- function(n, a = 1, b = 0, delta = 0.4) {
  nrm <- sqrt(a^2 + a * b + b^2 / 3)
  if (nrm == 0) stop("a and b cannot both be zero")
  T  <- runif(n)
  W  <- sqrt(T) * rnorm(n)
  cv <- (a * T + b * (T - T^2 / 2)) / nrm
  Vu <- cv / T * W + sqrt(pmax(1 - cv^2 / T, 0)) * rnorm(n)
  g  <- 0.5 + delta * sin(Vu)
  Y  <- rbinom(n, 1, g)
  f1 <- function(s) delta * exp(-0.5) * (a + b * (1 - s)) / nrm
  list(T = T, W = W, Y = Y, Vu = Vu, g = g, f1 = f1)
}
