# ogeom

**Nonlinear regression on sparsely observed Gaussian curves**

[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![R >= 3.5.0](https://img.shields.io/badge/R-%3E%3D%203.5.0-276DC3.svg)

`ogeom` implements Galerkin chaos estimators for binary regression on a Gaussian trajectory that is observed only at a few irregular times per subject, usually with measurement error. Beyond estimation, it answers a question that comes first: **given the times at which the curves are measured, which nonlinear effects of the trajectory can be learned at all, and how accurately?**

## Overview

In many longitudinal and biomedical studies the covariate of interest is a curve that is never seen in full: a biomarker profile, a growth trajectory or a physiological signal recorded at a handful of visits. For a binary response $Y$, the law of the data depends on the regression functional $g(X) = P(Y = 1 \mid X)$ only through the conditional expectation operator $K$, which maps a functional $h$ of the trajectory to $E\lbrace h(X) \mid O\rbrace$, where $O$ denotes the measurement times and the measurements of a subject. When $X$ is Gaussian, the Wiener–Itô chaos expansion gives

$$K^{*}K = E\lbrace \Gamma(B_D)\rbrace,$$

where $B_D$ is a random contraction determined by the design and the measurement error, and $\Gamma$ denotes second quantization. The information about effects of order $q$ is therefore the tensor moment $E(B_D^{\otimes q})$, while the information about linear effects is only its mean. Two designs with the same linear information can reveal very different nonlinear effects, and for Brownian motion and several other processes, effects of order $q$ are identified exactly when $q$ does not exceed the largest number of measurements per subject.

The package provides

- **exact conditional features** $\psi_\alpha(O) = E\lbrace b_\alpha(V) \mid O\rbrace$, Wick polynomials of the conditional means, with no plug-in approximation;
- **information operators** of a sampling design and their spectra, which show the effects a design can and cannot identify;
- **Galerkin estimators** of the chaos components, with an exact finite-sample risk, an unbiased risk criterion for the truncation level, and a bound-constrained version for laws estimated from data;
- **prediction** of $P(Y = 1 \mid O)$ and **first- and second-order kernels** for interpretation;
- **Fourier–Bessel estimators** for Brownian motion measured once at a uniformly distributed time;
- **data-generating functions** that reproduce all simulation designs of the accompanying paper.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("UfukBeyaztas/ogeom")
```

`devtools::install_github("UfukBeyaztas/ogeom")` works as well. The package `fdapace`, used only to estimate the law of the trajectory from real data, can be installed from CRAN:

```r
install.packages("fdapace")
```

## User manual

The reference manual documents every exported function: [**ogeom_0.1.0.pdf**](ogeom_0.1.0.pdf). After installation, the same help pages are available in R:

```r
help(package = "ogeom")
?og_fit
```

## Quick start

The interaction scenario of the paper uses a rank-four Gaussian process, between one and four noisy measurements per subject, and the regression function $g = 1/2 + 0.45\sin(V_1)\sin(V_2)$. This effect has no linear component, so linear methods cannot detect it, while its second-order component $0.45\,e^{-1}V_1V_2$ is identified because some subjects are measured at least twice.

```r
library(ogeom)
set.seed(1)
des   <- og_simulation_design("C", scenario = "interaction")
train <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
test  <- og_simulate(2000, des$process, des$design, des$sigma2, des$g, des$r)

# Galerkin chaos estimator of orders 1 and 2
fit <- og_fit(train, train$Y, des$process, des$sigma2, Q = 2, design = des$design)
fit                                   # truncation levels chosen by the risk criterion
round(fit$orders[[2]]$theta, 3)       # the V1*V2 coefficient estimates 0.45 * exp(-1) = 0.166

# prediction, and distance to the Bayes predictor P(Y = 1 | O)
p <- og_predict(fit, test)
sqrt(mean((p - test$Kg)^2))

# kernels of the fitted components
k <- og_kernels(fit)
image(k$grid, k$grid, k$beta2, xlab = "s", ylab = "t", main = "Second-order kernel")
```

## What a design can reveal

The information blocks of a design do not depend on the responses, so they can be computed before any data are collected. In the rotating example of the paper, $X(t) = Z_1\cos(2\pi t) + Z_2\sin(2\pi t)$ is measured once per subject. Design A places the time at $0$ or $1/4$, design B makes it uniform on $[0, 1]$. Both carry the same information about every linear functional of $X$:

```r
for (d in c("A", "B")) {
  des <- og_simulation_design("A", design = d)
  print(round(og_information_profile(des$process, des$design, des$sigma2, q = 2)$values, 3))
}
# design A: 0.5 0.5 0     the interaction Z1 Z2 is invisible
# design B: 0.5 0.25 0.25 every second-order effect is identified
```

## Simulation designs of the paper

| Experiment | What it shows | Data generation |
|---|---|---|
| A | equal linear information, different second-order information | `og_simulation_design("A", design = "A")` or `design = "B"`, then `og_simulate()` |
| B | the measurement threshold: an interaction identified only through subjects measured twice | `og_simulation_design("B", prob2 = 0.5)`, then `og_simulate()` |
| C | comparison with competing methods, interaction and logistic scenarios | `og_simulation_design("C", scenario = "interaction")` or `scenario = "usual"`, then `og_simulate()` |
| D | Brownian motion measured once: spectral and augmented Fourier–Bessel estimators | `og_simulate_bm(n, a, b)` |

Each call to `og_simulation_design()` returns the process, the design sampler, the error variance and the regression function, so a data set is generated with

```r
des <- og_simulation_design("B", prob2 = 0.5)
dat <- og_simulate(1000, des$process, des$design, des$sigma2, des$g, des$r)
```

`og_simulate()` draws the latent scores and the measurements from their exact joint Gaussian law, without discretizing the trajectory, and returns the Bayes predictor $P(Y = 1 \mid O)$ for evaluation.

## Brownian motion measured once

For Brownian motion measured once, without error, at a uniformly distributed time, the first-order information operator has kernel $-\log\max(s, u)$, eigenvalues $4/j_{0,k}^2$, where $j_{0,k}$ are the positive zeros of the Bessel function $J_0$, and Fourier–Bessel eigenfunctions. The spectral estimator attains the minimax rate $n^{-2b/(2b+3)}$ over Fourier–Bessel ellipsoids. When the first-order kernel does not vanish at the end of the interval, the augmented estimator restores that rate:

```r
set.seed(1)
d     <- og_simulate_bm(2000, a = 1, b = 0)      # constant kernel, nonzero at s = 1
J     <- og_bm_J(2000, b = 2)
fit_s <- og_bm_spectral(d$T, d$W, d$Y, J)
fit_a <- og_bm_augmented(d$T, d$W, d$Y, J)

s <- seq(0, 1, length.out = 201)
plot(s, d$f1(s), type = "l", lwd = 2, ylim = c(-0.1, 0.4), ylab = "first-order kernel")
lines(s, og_bm_kernel(s, fit_s), col = 2)
lines(s, og_bm_kernel(s, fit_a), col = 4)
legend("bottomleft", c("true", "spectral", "augmented"), col = c(1, 2, 4), lwd = c(2, 1, 1), bty = "n")
```

## Real data

With real data, the law of the trajectory is estimated from the pooled measurements, without the responses, and then treated as known. Here `Lt` and `Ly` are lists with the measurement times and the measurements of each subject, and `Y` holds the binary responses:

```r
library(fdapace)
fp   <- FPCA(Ly, Lt, list(dataType = "Sparse", error = TRUE, FVEthreshold = 0.99))
law  <- og_process_from_fpca(fp, sigma2_min = 0.01)
dat  <- og_center_data(list(t = Lt, w = Ly), law$mu_fun)
fit  <- og_fit_constrained(dat, Y, law$process, law$sigma2, Q = 2, gram = "conditional")
prob <- og_predict(fit, dat)
```

The chaos expansion is taken with respect to the Gaussian law. A skewed trajectory can make a linear effect look nonlinear, so it is worth transforming skewed biomarkers, for example by a Box–Cox transformation chosen without the responses, before fitting.

## Main functions

| Task | Functions |
|---|---|
| Law of the trajectory | `og_process_kl()`, `og_process_bm()`, `og_process_from_fpca()`, `og_center_data()`, `og_fourier_phi()`, `og_bm_dir_kl()`, `og_bm_dir_const()` |
| Conditional moments and features | `og_moments_all()`, `og_features()`, `og_multi_indices()`, `og_hermite_features()` |
| Information operators | `og_design_betas()`, `og_gram()`, `og_information_profile()` |
| Estimation and prediction | `og_fit()`, `og_fit_constrained()`, `og_predict()`, `og_component()`, `og_kernels()` |
| Brownian motion | `og_bm_spectral()`, `og_bm_augmented()`, `og_bm_kernel()`, `og_bm_J()`, `og_fb_eigen()`, `og_fb_phi()`, `og_bessel_j0_zeros()`, `og_bm_dir_fb()` |
| Simulation | `og_simulation_design()`, `og_simulate()`, `og_simulate_bm()`, `og_scenario()`, `og_design_sampler()`, `og_true_coefficients()`, `og_oracle_Kg()` |

## License

GPL (>= 3)

## Contact

Ufuk Beyaztas, Department of Statistics, Marmara University
[ufuk.beyaztas@marmara.edu.tr](mailto:ufuk.beyaztas@marmara.edu.tr)

Bug reports and suggestions are welcome through the [issue tracker](https://github.com/UfukBeyaztas/ogeom/issues).
