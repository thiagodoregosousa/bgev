# Backlog

## Estimator — done this version

`bgev_mle` now uses a multistart Nelder-Mead local search seeded by the
quantile start (`bgev_start_using_quantiles`) plus perturbed starts from a loose
data-driven box (`bgev_start_box`). The search runs on a reparametrised scale
(`log(sigma)`, `log(delta)`): estimation is restricted to `delta > 0` (bimodality
needs it, and `delta < 0` gives an unbounded likelihood at `x = mu`), and the
1e100 penalty only guards the data-dependent support; default `maxit` raised. It
returns point estimates, standard errors (`se`, from the inverse
observed-information Hessian, reported only for an admissible optimum) plus
diagnostics (convergence, multistart agreement, `admissible` PD-Hessian gate,
support-boundary check, `optimum`) and warns near the boundary and on tied data.
`likelihood = "grouped_likelihood"` added for discrete/rounded data.
`bgev_profile_likelihood` added. Dropped the grid-search / LHS region
(`bgev_start_region`) as redundant.

Validated by Monte Carlo (`benchmarks/mc_study.R`): regular regime (xi >= -0.2)
has ~0 bias, RMSE ~ 1/sqrt(n), near-nominal Wald coverage; the xi <~ -0.3
boundary is non-regular (Wald coverage collapses, flagged by `admissible`).
Methodology + results written up in `vignettes/bgev-estimation.Rmd`.

## Estimator — next version

- DE (DEoptim) as a *fallback* engine when the multistart disagrees or the
  quantile start fails (re-add DEoptim to `Imports`) — low priority, MC shows
  93-100% convergence without it.

Confidence intervals beyond the Hessian `se` are intentionally out of scope;
users can compute them with other packages if needed.

## Estimator — further validation (active research)

Validate against arXiv:2109.12738 Prop 3.8 (tail behaviour) and eq. 3.6
(quantile).


## CRAN resubmission

Version bump, `NEWS.md`, `cran-comments.md`, full `R CMD check --as-cran`,
submit.
