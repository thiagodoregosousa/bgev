# Consistency checks for continuous distribution implementations

Tests whether a set of functions implementing a continuous distribution
(density, distribution, quantile, and random generation) satisfy basic
probabilistic consistency conditions under the standard R naming
convention (`d*`, `p*`, `q*`, `r*`).

## Usage

``` r
distCheck(
  fun = "norm",
  n = 1000,
  robust = FALSE,
  subdivisions = 1500,
  support.lower = -Inf,
  support.upper = Inf,
  var.exists = TRUE,
  print.result = TRUE,
  ...
)
```

## Arguments

- fun:

  Character string giving the name of the distribution (e.g., `"norm"`,
  `"gev"`, `"exp"`).

- n:

  Sample size used when generating random values via the corresponding
  `r*` function.

- robust:

  Logical; if `TRUE`, mean and variance are computed using robust
  estimators when applicable.

- subdivisions:

  Number of subdivisions used for numerical integration when evaluating
  the density function.

- support.lower:

  Lower bound of the support of the distribution.

- support.upper:

  Upper bound of the support of the distribution.

- var.exists:

  Logical; indicates whether the variance of the distribution exists
  (useful for GEV, bimodal GEV, stable distributions, etc.).

- print.result:

  Logical; if `TRUE`, a summary of the test results is printed.

- ...:

  Additional parameters passed to the distribution functions.

## Value

A list containing the computed values, theoretical expectations, and
diagnostic information for each test.

## Details

This function is an adaptation of `fBasics::distCheck`, extended to
allow for custom distribution support and to return all test results in
a structured object for further inspection.

The following consistency checks are performed:

- Density check:

  Tests whether the density integrates to one over the specified
  support. For distributions with restricted support (e.g., GEV or
  bimodal GEV), appropriate bounds should be supplied.

- Quantile–CDF check:

  Compares empirical quantiles obtained from random generation with
  those implied by the cumulative distribution function.

- Mean–variance check:

  Computes mean and variance both from numerical integration of the
  density and from simulated samples, and compares the two. This check
  is skipped or flagged when moments are not finite.

## See also

`distCheck`

## Author

Thiago do Rego Sousa

## Examples

``` r
if (FALSE) { # \dontrun{
distCheck("norm")
distCheck("gev", xi = 0.2, sigma = 1, mu = 0,
          support.lower = -5, support.upper = 10)
} # }
```
