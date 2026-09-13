# Maximum Likelihood Estimation for the BGEV distribution

Maximum Likelihood Estimation for the BGEV distribution

## Usage

``` r
bgev_mle(
  x,
  lower = c(-12, 0.01, -12, -0.99),
  upper = c(12, 12, 12, 12),
  control = DEoptim::DEoptim.control(itermax = 100, NP = 100, trace = FALSE),
  DEoptim_replicates = 5
)
```

## Arguments

- x:

  Numeric vector of observations.

- lower:

  Optional vector of lower bounds for the parameters (mu, sigma, xi,
  delta).

- upper:

  Optional vector of upper bounds for the parameters (mu, sigma, xi,
  delta).

- control:

  List of control parameters, as returned by
  [DEoptim.control](https://rdrr.io/pkg/DEoptim/man/DEoptim.control.html).

- DEoptim_replicates:

  Number of independent DEoptim runs; the run with the best (highest)
  log-likelihood is returned.

## Value

The [DEoptim](https://rdrr.io/pkg/DEoptim/man/DEoptim.html) result for
the best-performing replicate: a list with `optim$bestmem` (the
estimated `c(mu, sigma, xi, delta)`) and `optim$bestval` (the negative
log-likelihood at that estimate).

## Author

Thiago do Rego Sousa and Yasmin Lirio

## Examples

``` r
# \donttest{
set.seed(1)
x <- rbgev(n = 200, mu = 1, sigma = 1, xi = 1, delta = 1)
fit <- bgev_mle(x, control = DEoptim::DEoptim.control(itermax = 20, NP = 40, trace = FALSE))
fit$optim$bestmem
#>      par1      par2      par3      par4 
#> 0.8509898 1.5183835 1.7577140 1.4065316 
# }
```
