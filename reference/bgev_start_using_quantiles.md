# Starting values for BGEV distribution

Internal helper: data-driven starting values `c(mu, sigma, xi, delta)`
from quantile matching, used to seed
[bgev_mle](https://thiagodoregosousa.github.io/bgev/reference/bgev_mle.md).

## Usage

``` r
bgev_start_using_quantiles(x)
```

## Arguments

- x:

  Numeric vector of observations.

## Value

A length-4 numeric vector of starting values for (mu, sigma, xi, delta).

## Author

Thiago do Rego Sousa
