# bgev 0.3

## Estimation

* `bgev_mle()` was rewritten. It now uses a multistart Nelder-Mead search on a
  reparametrised scale (`log(sigma)`, `log(delta)`), seeded by quantile matching
  plus loose data-driven box starts. Estimation is restricted to `delta > 0`
  (bimodality requires it, and for `delta < 0` the likelihood is unbounded at
  `x = mu`); the distribution functions still accept the full `delta > -1`.
* The return value changed: instead of the old `DEoptim` object
  (`$optim$bestmem`), `bgev_mle()` now returns `$par`, `$se`, `$loglik` and
  diagnostics. **This is a breaking change.**
* Standard errors (`se`) are returned from the inverse observed-information
  Hessian, for an admissible (regular) optimum only; near the parameter-dependent
  support boundary they are not reliable and are returned as `NA`.
* New diagnostics on every fit: `convergence`, `agree`, `admissible` (a
  positive-definite-Hessian acceptance gate that rejects spurious optima),
  `boundary`, and `optimum`.
* `likelihood = "grouped_likelihood"` added for discrete or rounded data, using
  the interval likelihood `F(x + h/2) - F(x - h/2)`; tied data under the
  continuous density now raises a warning.
* New `bgev_profile_likelihood()` for the profile log-likelihood of a parameter.

## Documentation

* New vignette "Maximum Likelihood Estimation for the BGEV Distribution",
  including a Monte Carlo validation of the estimator.

## Internal

* `distCheck()` renamed to `dist_check()`.
* Imports: dropped `DEoptim` and `lhs`; added `numDeriv` and `graphics`.
