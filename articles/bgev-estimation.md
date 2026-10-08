# Maximum Likelihood Estimation for the BGEV Distribution

This vignette documents how
[`bgev_mle()`](https://thiagodoregosousa.github.io/bgev/reference/bgev_mle.md)
estimates the bimodal generalized extreme value (BGEV) distribution: the
likelihood and a singularity that forces a restriction on `delta`, the
starting-value strategy, a grouped likelihood for discrete data, the
diagnostics returned with every fit, and a Monte Carlo study validating
the estimator.

The parametrization follows the *revised* BGEV of Otiniano, Lisboa &
Ribeiro (2025), which adds a location parameter. With shape `xi`,
location `mu`, scale `sigma > 0` and bimodality parameter `delta > -1`,
the CDF is
$`F(x;\xi,\mu,\sigma,\delta) = F_{\mathrm{GEV}}(T(x);\xi,0,\sigma)`$
with the transform $`T(x) = (x-\mu)\,|x-\mu|^{\delta}`$ and derivative
$`T'(x) = (\delta+1)\,|x-\mu|^{\delta}`$.

## 1. The likelihood and its singularity at `x = mu`

The density is $`f(x) = f_{\mathrm{GEV}}(T(x);\xi,0,\sigma)\,T'(x)`$.
The factor $`T'(x) = (\delta+1)|x-\mu|^{\delta}`$ behaves very
differently with the sign of `delta`:

- for $`\delta > 0`$, $`T'(x) \to 0`$ as $`x \to \mu`$ (the bimodality
  *dip*);
- for $`\delta < 0`$, $`T'(x) \to +\infty`$ as $`x \to \mu`$ — **the
  density diverges at $`x = \mu`$.**

``` r

near_mu <- 10^-(1:4)
rbind(
  `delta<0` = sapply(near_mu, function(e) dbgev(2 + e, mu = 2, sigma = 1, xi = -0.3, delta = -0.3)),
  `delta>0` = sapply(near_mu, function(e) dbgev(2 + e, mu = 2, sigma = 1, xi = -0.3, delta =  0.5))
)
#>              [,1]       [,2]       [,3]        [,4]
#> delta<0 0.5358177 1.03675932 2.05034061 4.083283699
#> delta>0 0.1760839 0.05519845 0.01745022 0.005518193
```

When `delta < 0` the likelihood is therefore **unbounded**: driving `mu`
onto an observation sends the log-likelihood to $`+\infty`$. This is the
classic unbounded-likelihood phenomenon — the same mechanism as the
normal model with $`\sigma \to 0`$ on a single point (Pawitan, 2001,
§4.8). With continuous data it has probability zero, but with tied or
integer data the global maximizer parks `mu` on a data atom, a spurious
spike rather than a better fit.

**Consequence for estimation.**
[`bgev_mle()`](https://thiagodoregosousa.github.io/bgev/reference/bgev_mle.md)
restricts estimation to `delta > 0`. This matches the revised reference,
whose own estimation routine uses `delta >= 0`, and it loses nothing of
interest: bimodality occurs only for `delta > 0` (`delta = 0` is the
ordinary GEV). The *distribution* functions
([`dbgev()`](https://thiagodoregosousa.github.io/bgev/reference/bgev.md),
[`pbgev()`](https://thiagodoregosousa.github.io/bgev/reference/bgev.md),
…) still accept the full `delta > -1`; only the estimator is restricted.
Internally the search is reparametrized as
$`(\mu, \log\sigma, \log\delta)`$, so `sigma > 0` and `delta > 0` hold
by construction and no penalty cliff is needed for the parameter box.

## 2. Starting values: quantiles for the aim, a box for diversity

The optimizer is a **multistart local search** (Nelder–Mead, 1965). Two
complementary kinds of starting point feed it:

1.  **One quantile-matching start** (`bgev_start_using_quantiles`), the
    *aimed* guess: it finds parameters whose model quantiles match the
    empirical quantiles at $`p = (0.1, 0.3, 0.6, 0.9)`$ by solving the
    four equations with a Newton-type solver.
2.  **Several box starts** (`bgev_start_box`, `bgev_sample_starts`), the
    *diversity*: points drawn uniformly from a loose, data-driven box.
    They are deliberately scattered to probe for competing local maxima.

[`bgev_mle()`](https://thiagodoregosousa.github.io/bgev/reference/bgev_mle.md)
runs a full local optimization from every start and keeps the best
*admissible* result (Section 4). The aimed start usually wins; the box
starts guard against multimodality.

This local-multistart design is **cheaper than a global optimizer** such
as differential evolution (`DEoptim`). The accepted pattern is to find
the basin with cheap starts and then polish locally (the global-to-local
strategy; Nocedal & Wright, 2006; Nash, 2014): the revised reference
fits with plain Nelder–Mead, and the Monte Carlo study below reaches the
maximum in 93–100% of replications, so a full global search is
unnecessary for the regular regime. A global fallback is reserved for
cases where the starts disagree.

### Closed-form quantile estimators

The BGEV quantile function yields exact estimators at special
probabilities. Writing $`Q(p)`$ for the quantile function, at
$`p = e^{-1}`$ we have $`-\log p = 1`$, so the shape term
$`[(-\log p)^{-\xi}-1]/\xi`$ vanishes and

``` math
Q(e^{-1}) = \mu \qquad \text{for all } \sigma, \xi, \delta.
```

``` r

x <- rbgev(2e5, mu = 5, sigma = 3, xi = -0.4, delta = 0.3)
c(closed_form = unname(quantile(x, exp(-1))), truth = 5)
#> closed_form       truth 
#>    5.010103    5.000000
```

So `quantile(x, exp(-1))` is an exact, parameter-free estimator of `mu`.
For the Gumbel-type case `xi = 0`, closed forms exist for the other two
parameters as well. With $`q_1 = Q(e^{-e^{2}})`$ and
$`q_2 = Q(e^{-e^{1}})`$,

``` math
\delta = \frac{1}{\log_2(q_1/q_2)} - 1, \qquad \sigma = (-q_2)^{\delta+1}.
```

These are documented for completeness. They are **not** used as the
estimator: the general `xi` case has no such closed form, and
experiments showed that seeding `mu` with $`Q(e^{-1})`$ did not improve
convergence or recovery over the joint quantile solve — starts were not
the bottleneck.

## 3. Discrete data: the grouped likelihood

BGEV is continuous, so a fit to rounded or integer data via the density
is a model mismatch (and, were `delta < 0` allowed, exactly where the
atom spike bites). The correct model records each value `x` to
resolution `h` as the interval probability

``` math
P\big(X \in [x - h/2,\ x + h/2]\big) = F(x + h/2) - F(x - h/2),
```

which is bounded by one and cannot blow up (Pawitan, 2001, §4.8). Select
it with the `likelihood` argument:

``` r

xi_int <- sample(10:40, 120, replace = TRUE)       # integer data
fit_grp <- bgev_mle(xi_int, likelihood = "grouped_likelihood", h = 1)
round(fit_grp$par, 3)
#>     mu  sigma     xi  delta 
#> 21.269 14.542 -0.168  0.220
```

On continuous data the grouped likelihood with a small `h` reproduces
the continuous MLE; on discrete data it is the appropriate choice and
suppresses the tied-data warning.

## 4. Diagnostics returned with every fit

``` r

x   <- rbgev(300, mu = 0, sigma = 1, xi = 0.2, delta = 1)
fit <- bgev_mle(x)
round(fit$par, 3)
#>     mu  sigma     xi  delta 
#> -0.004  0.918  0.315  1.232
c(loglik = fit$loglik, convergence = fit$convergence,
  admissible = fit$admissible, agree = fit$agree)
#>      loglik convergence  admissible       agree 
#>   -340.5378      0.0000      1.0000      0.0000
```

- **`convergence`** — the `optim` code (0 = success).
- **`agree`** — whether several independent starts reached the same
  maximum (a practical “is it global?” signal).
- **`admissible`** — whether the returned optimum is a *regular*
  interior maximum: a positive-definite Hessian with eigenvalue ratio
  above `1e-6`. The multistart keeps the best admissible optimum and
  rejects spurious / singular ones; `admissible = FALSE` warns that
  inference is not trustworthy there.
- **`boundary`** — proximity of the data to the parameter-dependent
  support endpoint, where the usual regularity conditions fail.
- **`optimum`** — gradient norm, Hessian positive-definiteness and
  eigenvalue ratio at the estimate.

[`bgev_profile_likelihood()`](https://thiagodoregosousa.github.io/bgev/reference/bgev_profile_likelihood.md)
complements these with a profile curve for any parameter — the check
Otiniano et al. use to confirm an optimum is global.

## 5. Monte Carlo validation

A study of 13,500 fits (45 cells, 300 replications) simulated from known
parameters and refit, with `mu = 0`, `sigma = 1`, `xi` in {−0.4, −0.2,
0, 0.2, 0.4}, `delta` in {0.25, 1, 3} and `n` in {100, 250, 500}. The
script is `benchmarks/mc_study.R`. Convergence was 93–100% in every cell
(maximum failure rate 6%). The results below summarize the two
decision-relevant findings; they are static, taken from that run.

The estimator is **consistent and well-calibrated in the regular
regime** (`xi >= -0.2`): bias → 0, RMSE falls like `1/sqrt(n)`, and Wald
coverage from the observed-information Hessian is near the nominal 0.95.

|   xi | n=100 | n=250 | n=500 |
|-----:|------:|------:|------:|
| -0.4 |  0.84 |  0.39 |  0.08 |
| -0.2 |  0.98 |  0.97 |  0.97 |
|  0.0 |  0.92 |  0.94 |  0.95 |
|  0.2 |  0.91 |  0.94 |  0.93 |
|  0.4 |  0.91 |  0.93 |  0.94 |

Wald 95% coverage of xi (delta = 1). Target 0.95. {.table}

At the boundary (`xi = -0.4`) the model is **non-regular**: the
parameter-dependent support makes the MLE non-normal, so Wald coverage
collapses *and worsens with n* (0.84 → 0.39 → 0.08) and RMSE does not
shrink. The package flags this itself through the admissibility gate,
which rejects most such fits:

|   xi | delta=0.25 | delta=1 | delta=3 |
|-----:|-----------:|--------:|--------:|
| -0.4 |       0.46 |    0.12 |    0.10 |
| -0.2 |       0.98 |    0.92 |    0.89 |
|  0.0 |       1.00 |    1.00 |    1.00 |
|  0.2 |       1.00 |    1.00 |    1.00 |
|  0.4 |       1.00 |    1.00 |    1.00 |

Admissible rate (positive-definite-Hessian gate), n = 500. {.table}

Where Wald coverage holds the gate accepts ~100% of fits; where it fails
it rejects 88–90%. A user who checks `fit$admissible` is steered away
from exactly the fits whose Hessian-based standard errors cannot be
trusted. For those boundary cells, parametric-bootstrap or
profile-likelihood intervals are the appropriate alternative to Wald
intervals.

## References

- Otiniano, C. E. G., Lisboa, M. N. S., & Ribeiro, T. K. A. (2025). *A
  Revised Bimodal Generalized Extreme Value Distribution: Theory and
  Climate Data Application.* Entropy, 27(7), 749.
  [doi:10.3390/e27070749](https://doi.org/10.3390/e27070749)
- Otiniano, C. E. G., et al. (2023). *A bimodal model for extremes
  data.* Environmental and Ecological Statistics.
  [doi:10.1007/s10651-023-00566-7](https://doi.org/10.1007/s10651-023-00566-7)
- Pawitan, Y. (2001). *In All Likelihood: Statistical Modelling and
  Inference Using Likelihood.* Oxford University Press. (§4.8, unbounded
  likelihood and the finite-precision / grouped likelihood.)
- Nelder, J. A., & Mead, R. (1965). *A simplex method for function
  minimization.* The Computer Journal, 7(4), 308–313.
- Nash, J. C. (2014). *Nonlinear Parameter Optimization Using R Tools.*
  Wiley.
- Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd
  ed.). Springer.
