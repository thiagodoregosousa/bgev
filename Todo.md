# Backlog

## Estimator work (active research)

Finish the grid-search + bounded-region MLE in `R/bgev_start_params.R` /
`R/bgev_estimation.R`. Validate against arXiv:2109.12738 Prop 3.8 (tail
behaviour) and eq. 3.6 (quantile). Extend
`tests/testthat/test_bgev_estimation.R` / `test_bgev_start_quantile.R`
accordingly, then re-run `benchmarks/monte_carlo_study.R` once stable.

## Data files

`data/densidade_ar_max.xlsx` / `data/umidade_min_ex.csv` need to move to
`inst/extdata` and be documented (or be dropped if unused) — their
current placement in `data/` will trip `R CMD check`, which expects
`.rda`/`.RData` there.

## Confidence intervals / moments

Flesh out or drop `to_be_implemented/bgev_conf_intervals_and_moments.R`.

## CRAN resubmission

Version bump, `NEWS.md`, `cran-comments.md`, full
`R CMD check --as-cran`, submit.

## qqplot for bgev

No dedicated feature — at most, add a small qqplot demo to `qbgev`’s
`@examples` (generic Q-Q machinery already works against any `q*`
function, nothing bgev-specific to build).

## pkgdown site

Minimal `_pkgdown.yml` + GitHub Actions workflow to build and deploy to
`gh-pages`. Once scaffolded, group the reference index (distribution vs.
estimation vs. diagnostics) and theme it minimally.
