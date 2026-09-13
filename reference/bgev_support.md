# Compute the support of the BGEV distribution

Returns the lower and upper limits of the support of the BGEV
distribution

## Usage

``` r
bgev_support(mu = 1, sigma = 1, xi = 0.3, delta = 2)
```

## Arguments

- mu:

  location parameter

- sigma:

  scale parameter (sigma \> 0)

- xi:

  shape parameter in R

- delta:

  shape parameter (delta \> -1)

## Value

A vector of length 2 with the lower and upper limits of the support

## Details

It returns values with `-Inf` or `Inf` when the support is unbounded.
When the shape parameter `xi` is different from zero, the support is
truncated either at the left or at the right side of the real.
Considering the support is particularly useful to estimating momoments
and to compute the likelihood function.

## Author

Thiago do Rego Sousa
