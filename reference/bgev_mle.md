# Maximum Likelihood Estimation for the BGEV distribution

Fits the BGEV distribution by maximum likelihood with a multistart local
search (Nelder-Mead). The primary start comes from
[bgev_start_using_quantiles](https://thiagodoregosousa.github.io/bgev/reference/bgev_start_using_quantiles.md);
additional starts are drawn inside a loose data-driven box (see
`bgev_start_box`). The run with the highest log-likelihood is returned.
Multistart guards against the multiple local maxima that are known to
occur in bimodal likelihoods.

## Usage

``` r
bgev_mle(
  x,
  likelihood = c("continuous_density", "grouped_likelihood"),
  h = 1,
  n_starts = 10,
  k = 3,
  control = list(maxit = 2000),
  ...
)
```

## Arguments

- x:

  Numeric vector of observations.

- likelihood:

  Which likelihood to maximise: `"continuous_density"` (the default,
  using the density) or `"grouped_likelihood"` for discrete or rounded
  data, which uses the interval probability `F(x + h/2) - F(x - h/2)` of
  each observation (see `h`).

- h:

  Rounding resolution for `"grouped_likelihood"` (default 1, for integer
  data). Ignored when `likelihood = "continuous_density"`.

- n_starts:

  Number of starting points for the multistart search (the quantile
  start plus `n_starts - 1` perturbed starts).

- k:

  Half-width multiplier for the data-driven box of perturbed starts.

- control:

  List of control parameters passed to
  [optim](https://rdrr.io/r/stats/optim.html) (defaults to a higher
  `maxit` than `optim`'s own, which the bounded search otherwise hits on
  this likelihood).

- ...:

  Additional arguments passed to
  [optim](https://rdrr.io/r/stats/optim.html).

## Value

A list with `par` (named estimate `c(mu, sigma, xi, delta)`), `se`
(standard errors from the inverse observed-information Hessian, `NA`
unless `admissible`), `loglik` (maximised log-likelihood, positive, on
the chosen scale), `likelihood` (which likelihood was used),
`convergence` (`optim` code, 0 = success), `start` (the quantile start),
`n_starts`, `loglik_starts`, `agree` (TRUE when several starts reach the
selected maximum), `admissible` (TRUE when the returned optimum has a
positive-definite Hessian), `boundary` (support-boundary diagnostic),
and `optimum` (gradient norm, Hessian positive-definiteness, eigenvalue
ratio and `se`).

## Details

The fit carries diagnostics (convergence, agreement across starts, and a
support-boundary check) because the BGEV support depends on the
parameters, so the usual regularity conditions can fail near the
boundary. Standard errors from the inverse observed-information Hessian
are returned, but only for an admissible (regular) optimum; near the
boundary they are not reliable and are returned as `NA`.

Although the BGEV distribution is defined for `delta > -1`, estimation
is restricted to `delta > 0`. This is deliberate: bimodality – the
purpose of the model – occurs only for `delta > 0`, and for `delta < 0`
the density is unbounded at `x = mu` (the transform derivative
`(delta + 1)|x - mu|^delta` diverges), giving a singular likelihood
whose global maximum is a spurious spike on a data point. The revised
BGEV reference estimates on `delta >= 0` for the same reason.

The search is run on a reparametrised scale – `log(sigma)` and
`log(delta)` – so that `sigma > 0` and `delta > 0` hold automatically.
The only remaining penalty guards the data-dependent support
(observations beyond the fitted endpoint), which is not a box constraint
and cannot be transformed away. Among the multistart results, the fit
returned is the highest-likelihood one whose Hessian at the optimum is
positive-definite (a genuine interior maximum, rejecting spurious
spikes); `admissible` is `FALSE` if none qualified. Estimates are
reported on the natural scale.

## Author

Thiago do Rego Sousa and Yasmin Lirio

## Examples

``` r
# \donttest{
set.seed(1)
x <- rbgev(n = 200, mu = 1, sigma = 1, xi = 1, delta = 1)
fit <- bgev_mle(x)
#> Warning: No start produced a regular (positive-definite, non-singular) Hessian; the returned estimate may be a boundary or spurious optimum -- see $optimum and $boundary.
fit$par
#>        mu     sigma        xi     delta 
#> 0.9970241 0.8855053 0.9229374 1.0101225 
# }
```
