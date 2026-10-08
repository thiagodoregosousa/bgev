# Profile log-likelihood for a BGEV parameter

Computes the profile log-likelihood of one BGEV parameter over a grid,
maximising over the remaining three at each grid point (Nelder-Mead).
This is the diagnostic used in Otiniano et al. (2023) to check whether a
fitted optimum is global and the parameter is well identified.

## Usage

``` r
bgev_profile_likelihood(x, par, which, span = 0.5, n = 41, plot = TRUE)
```

## Arguments

- x:

  Numeric vector of observations.

- par:

  Vector `c(mu, sigma, xi, delta)`, e.g. `bgev_mle(x)$par`.

- which:

  Index (1-4) or name of the parameter to profile.

- span:

  Half-width of the grid, as a fraction of the parameter value.

- n:

  Number of grid points.

- plot:

  Logical; if `TRUE`, plot the profile curve.

## Value

A data frame with the grid values of the profiled parameter and the
profile log-likelihood.

## Author

Thiago do Rego Sousa
