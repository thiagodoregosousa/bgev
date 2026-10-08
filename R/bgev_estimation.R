#' Log-likelihood function for the BGEV distribution
#'
#' @param x Numeric vector of observations.
#' @param pars Vector of parameters (mu, sigma, xi, delta). See \link{bgev}.
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#'
#' @return The log-likelihood value.
#' @importFrom stats setNames
#' @export 
bgev_log_likelihood <- function(x, pars) {
  
  stopifnot(is.numeric(pars))
  stopifnot(length(pars) == 4)
  
  pars <- setNames(pars, c("mu", "sigma", "xi", "delta"))
  
  mu    <- pars["mu"]
  sigma <- pars["sigma"]
  xi    <- pars["xi"]
  delta <- pars["delta"]
  
  if (!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma > 0 and delta > -1")
  
  sum(log(dbgev(x, mu, sigma, xi, delta)))
}


# Negative log-likelihood used as the objective for the optimizers. It wraps the
# positive bgev_log_likelihood and returns a large finite penalty on invalid or
# non-finite points, so the search never leaves the valid parameter domain.
bgev_negative_log_likelihood <- function(x, pars) {
  if (!bgev_valid_params(pars[1], pars[2], pars[3], pars[4]))
    return(1e+100)
  val <- -bgev_log_likelihood(x, pars)
  if (!is.finite(val)) 1e+100 else val
}


# Grouped / interval negative log-likelihood for discrete (rounded, tied) data:
# a value recorded to resolution h contributes the interval probability
# F(x + h/2) - F(x - h/2) instead of the density f(x). An integrable singularity
# has finite mass, so this cannot blow up -- it is the correct model for discrete
# data and removes the delta < 0 density spike at x = mu. Ref: Pawitan, sec 4.8.
bgev_negative_log_likelihood_grouped <- function(x, pars, h) {
  if (!bgev_valid_params(pars[1], pars[2], pars[3], pars[4]))
    return(1e+100)
  p_hi <- pbgev(x + h / 2, pars[1], pars[2], pars[3], pars[4])
  p_lo <- pbgev(x - h / 2, pars[1], pars[2], pars[3], pars[4])
  val <- -sum(log(p_hi - p_lo))
  if (!is.finite(val)) 1e+100 else val
}


#' Maximum Likelihood Estimation for the BGEV distribution
#'
#' Fits the BGEV distribution by maximum likelihood with a multistart local
#' search (Nelder-Mead). The primary start comes from
#' \link{bgev_start_using_quantiles}; additional starts are drawn inside a loose
#' data-driven box (see \code{bgev_start_box}). The run with the highest
#' log-likelihood is returned. Multistart guards against the multiple local
#' maxima that are known to occur in bimodal likelihoods.
#'
#' The fit carries diagnostics (convergence, agreement across starts, and a
#' support-boundary check) because the BGEV support depends on the parameters,
#' so the usual regularity conditions can fail near the boundary. Standard errors
#' from the inverse observed-information Hessian are returned, but only for an
#' admissible (regular) optimum; near the boundary they are not reliable and are
#' returned as \code{NA}.
#'
#' @param x Numeric vector of observations.
#' @param likelihood Which likelihood to maximise: \code{"continuous_density"}
#'   (the default, using the density) or \code{"grouped_likelihood"} for discrete
#'   or rounded data, which uses the interval probability
#'   \code{F(x + h/2) - F(x - h/2)} of each observation (see \code{h}).
#' @param h Rounding resolution for \code{"grouped_likelihood"} (default 1, for
#'   integer data). Ignored when \code{likelihood = "continuous_density"}.
#' @param n_starts Number of starting points for the multistart search
#'   (the quantile start plus \code{n_starts - 1} perturbed starts).
#' @param k Half-width multiplier for the data-driven box of perturbed starts.
#' @param control List of control parameters passed to \link[stats]{optim}
#'   (defaults to a higher \code{maxit} than \code{optim}'s own, which the
#'   bounded search otherwise hits on this likelihood).
#' @param ... Additional arguments passed to \link[stats]{optim}.
#'
#' @details
#' Although the BGEV distribution is defined for \code{delta > -1}, estimation is
#' restricted to \code{delta > 0}. This is deliberate: bimodality -- the purpose
#' of the model -- occurs only for \code{delta > 0}, and for \code{delta < 0} the
#' density is unbounded at \code{x = mu} (the transform derivative
#' \code{(delta + 1)|x - mu|^delta} diverges), giving a singular likelihood whose
#' global maximum is a spurious spike on a data point. The revised BGEV reference
#' estimates on \code{delta >= 0} for the same reason.
#'
#' The search is run on a reparametrised scale -- \code{log(sigma)} and
#' \code{log(delta)} -- so that \code{sigma > 0} and \code{delta > 0} hold
#' automatically. The only remaining penalty guards the data-dependent support
#' (observations beyond the fitted endpoint), which is not a box constraint and
#' cannot be transformed away. Among the multistart results, the fit returned is
#' the highest-likelihood one whose Hessian at the optimum is positive-definite
#' (a genuine interior maximum, rejecting spurious spikes); \code{admissible} is
#' \code{FALSE} if none qualified. Estimates are reported on the natural scale.
#'
#' @return A list with \code{par} (named estimate \code{c(mu, sigma, xi, delta)}),
#'   \code{se} (standard errors from the inverse observed-information Hessian,
#'   \code{NA} unless \code{admissible}), \code{loglik} (maximised log-likelihood,
#'   positive, on the chosen scale), \code{likelihood} (which likelihood was
#'   used), \code{convergence} (\code{optim} code, 0 = success), \code{start}
#'   (the quantile start), \code{n_starts}, \code{loglik_starts}, \code{agree}
#'   (TRUE when several starts reach the selected maximum), \code{admissible}
#'   (TRUE when the returned optimum has a positive-definite Hessian),
#'   \code{boundary} (support-boundary diagnostic), and \code{optimum} (gradient
#'   norm, Hessian positive-definiteness, eigenvalue ratio and \code{se}).
#'
#' @author Thiago do Rego Sousa and Yasmin Lirio
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- rbgev(n = 200, mu = 1, sigma = 1, xi = 1, delta = 1)
#' fit <- bgev_mle(x)
#' fit$par
#' }
#' @importFrom stats optim runif
#' @export
bgev_mle <- function(x, likelihood = c("continuous_density", "grouped_likelihood"),
                     h = 1, n_starts = 10, k = 3, control = list(maxit = 2000), ...) {
  likelihood <- match.arg(likelihood)
  if (is.null(x) || !is.numeric(x) || anyNA(x))
    stop("`x` must be a non-null numeric vector with no missing values.", call. = FALSE)

  # Objective (natural scale). For discrete/rounded data use the grouped
  # (interval) likelihood; otherwise the continuous density. Tied values under
  # the continuous density are a modelling mismatch, so warn and point to the
  # grouped option.
  if (likelihood == "grouped_likelihood") {
    nll <- function(pars) bgev_negative_log_likelihood_grouped(x, pars, h)
  } else {
    nll <- function(pars) bgev_negative_log_likelihood(x, pars)
    if (anyDuplicated(x))
      warning("`x` has tied (repeated) values; BGEV is continuous. For rounded/",
              "discrete data use likelihood = \"grouped_likelihood\".", call. = FALSE)
  }

  # Search on the reparametrised scale: theta = (mu, log(sigma), xi, log(delta)).
  # sigma > 0 and delta > 0 hold by construction, so Nelder-Mead never meets the
  # parameter-box cliff and the delta < 0 singularity is unreachable -- only the
  # data-dependent support penalty remains.
  neg_loglik_theta <- function(theta) nll(bgev_from_theta(theta))

  # Quantile start lands in the basin; the loose box feeds perturbed starts.
  # Every start needs delta > 0 before the log(delta) transform.
  delta_floor <- 1e-2
  theta0 <- bgev_start_using_quantiles(x)
  theta0[4] <- max(theta0[4], delta_floor)
  box <- bgev_start_box(x, k = k)
  starts <- rbind(theta0, bgev_sample_starts(n_starts - 1, box))
  starts[, 4] <- pmax(starts[, 4], delta_floor)

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    stats::optim(par = bgev_to_theta(starts[i, ]), fn = neg_loglik_theta,
                 method = "Nelder-Mead", control = control, ...)
  })

  values <- vapply(fits, function(f) f$value, numeric(1))  # negative log-lik
  loglik_starts <- -values

  # Acceptance gate: walk candidates from best to worst log-likelihood and keep
  # the first whose Hessian at the optimum is positive-definite -- a genuine
  # interior maximum. This rejects spurious spikes (non-PD / singular Hessian).
  # If none qualifies, fall back to the best and flag it with a warning.
  ord <- order(values)
  selected <- NULL; diag_sel <- NULL
  for (i in ord) {
    diag_i <- bgev_optimum_diagnostics(nll, bgev_from_theta(fits[[i]]$par))
    if (bgev_is_regular_optimum(diag_i)) { selected <- fits[[i]]; diag_sel <- diag_i; break }
  }
  admissible <- !is.null(selected)
  if (!admissible) {
    selected <- fits[[ord[1]]]
    diag_sel <- bgev_optimum_diagnostics(nll, bgev_from_theta(selected$par))
    warning("No start produced a regular (positive-definite, non-singular) ",
            "Hessian; the returned estimate may be a boundary or spurious ",
            "optimum -- see $optimum and $boundary.", call. = FALSE)
  }

  par <- stats::setNames(bgev_from_theta(selected$par), c("mu", "sigma", "xi", "delta"))
  selected_loglik <- -selected$value
  # Several starts reaching the selected maximum is the practical signal that it
  # is global; disagreement flags a multimodal likelihood.
  agree <- sum(abs(loglik_starts - selected_loglik) < 1e-3) > 1

  boundary <- bgev_check_boundary(x, par)
  if (boundary$near)
    warning("Estimate is near the support boundary; the standard errors from ",
            "the Hessian are not reliable here (regularity conditions fail).",
            call. = FALSE)

  # Standard errors from the inverse observed-information Hessian. Only reported
  # for an admissible (regular) optimum; near the parameter-dependent support
  # boundary the Hessian is not trustworthy, so `se` is NA there.
  pnames <- c("mu", "sigma", "xi", "delta")
  se <- if (admissible) diag_sel$se else rep(NA_real_, 4)
  se <- stats::setNames(se, pnames)

  list(par = par, se = se, loglik = selected_loglik, likelihood = likelihood,
       convergence = selected$convergence,
       start = stats::setNames(theta0, pnames),
       n_starts = nrow(starts), loglik_starts = loglik_starts,
       agree = agree, admissible = admissible, boundary = boundary,
       optimum = diag_sel)
}


# Reparametrisation helpers. theta = (mu, log(sigma), xi, log(delta)). Estimation
# restricts delta > 0 (see bgev_mle details), so log(delta) is well defined and
# the delta < 0 density singularity at x = mu is unreachable by construction.
bgev_to_theta <- function(par)   c(par[1], log(par[2]), par[3], log(par[4]))
bgev_from_theta <- function(theta) c(theta[1], exp(theta[2]), theta[3], exp(theta[4]))


# Optimum diagnostics (KKT-style): gradient norm and Hessian of the negative
# log-likelihood at the estimate, on the natural scale. A positive-definite
# Hessian with a healthy eigenvalue ratio means a genuine, well-identified
# maximum; a near-singular Hessian (tiny ratio) is the flat ridge seen when the
# model is weakly identified (e.g. delta ~ 0, reducing to the GEV). Evaluated
# near the support boundary the finite differences can be unreliable, so the
# whole thing is wrapped and returns NA on failure.
# A regular (trustworthy) interior optimum has a positive-definite Hessian that
# is not effectively singular; an eigenvalue ratio below ~1e-6 signals a flat
# direction (optimx uses the same threshold). Used as the acceptance gate.
bgev_is_regular_optimum <- function(d, eig_ratio_tol = 1e-6) {
  isTRUE(d$hessian_pd) && is.finite(d$eig_ratio) && d$eig_ratio > eig_ratio_tol
}


bgev_optimum_diagnostics <- function(nll, par) {
  na <- list(grad_norm = NA_real_, hessian_pd = NA, eig_min = NA_real_,
             eig_max = NA_real_, eig_ratio = NA_real_, se = rep(NA_real_, length(par)))
  if (!requireNamespace("numDeriv", quietly = TRUE)) return(na)
  tryCatch({
    g <- numDeriv::grad(nll, par)
    H <- numDeriv::hessian(nll, par)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    # Standard errors from the inverse observed-information Hessian (H is the
    # Hessian of the NEGATIVE log-likelihood, so H itself is the observed
    # information). Valid only where the optimum is regular; a negative variance
    # is returned as NA.
    se <- tryCatch({
      v <- diag(solve(H)); v[v < 0] <- NA_real_; sqrt(v)
    }, error = function(e) rep(NA_real_, length(par)))
    list(grad_norm = sqrt(sum(g^2)),
         hessian_pd = all(eig > 0),
         eig_min = min(eig), eig_max = max(eig),
         eig_ratio = min(eig) / max(eig), se = se)
  }, error = function(e) na)
}


# Support-boundary diagnostic. The BGEV support endpoint depends on the
# parameters, so when the data reach close to the fitted endpoint the regularity
# conditions fail and likelihood-based standard errors are untrustworthy. We
# flag proximity relative to the data scale (a heuristic margin).
bgev_check_boundary <- function(x, par) {
  support <- bgev_support(par[1], par[2], par[3], par[4])
  scale <- stats::mad(x)
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(x)

  side <- NA_character_; endpoint <- NA_real_; margin <- Inf
  if (par[3] > 0) {          # finite lower endpoint (Frechet-type)
    side <- "lower"; endpoint <- support[1]; margin <- min(x) - endpoint
  } else if (par[3] < 0) {   # finite upper endpoint (Weibull-type)
    side <- "upper"; endpoint <- support[2]; margin <- endpoint - max(x)
  }
  near <- is.finite(margin) && margin < 0.05 * scale
  list(near = near, side = side, endpoint = endpoint, margin = margin)
}
























#' Starting values for BGEV distribution
#'
#' Internal helper: data-driven starting values \code{c(mu, sigma, xi, delta)}
#' from quantile matching, used to seed \link{bgev_mle}.
#'
#' @param x Numeric vector of observations.
#' @return A length-4 numeric vector of starting values for (mu, sigma, xi, delta).
#' @author Thiago do Rego Sousa
#' @keywords internal
#' @importFrom stats mad median quantile sd
bgev_start_using_quantiles = function(x){
  # Note: mu has an exact closed form, Q(e^-1) = mu (revised BGEV quantile).
  # Seeding mu with it and solving only 3 equations was tested against the full
  # 4-equation solve on the small-n and boundary MC cells and gave no improvement
  # in convergence/failure rate and slightly worse mu RMSE, so the joint 4-eq
  # solve (which the Monte Carlo study validated) is kept.
  p <- c(0.1, 0.3, 0.6, 0.9)
  q_emp <- quantile(x, p)

  quantile_error <- function(theta){
    mu  <- theta[1]
    sigma <- theta[2]
    xi  <- theta[3]
    delta <- theta[4]

    # The solver can probe invalid parameters (sigma <= 0, delta <= -1) where
    # qbgev is undefined; return a large residual there to steer it back instead
    # of letting qbgev error out.
    if (!bgev_valid_params(mu, sigma, xi, delta))
      return(rep(1e6, length(p)))

    q_model <- qbgev(p, mu, sigma, xi, delta)
    q_model - q_emp
  }

  start <- c(median(x), mad(x), 0.1, 1) # starting values for mu, sd, xi, delta
  sol <- nleqslv::nleqslv(start, quantile_error)$x # solve nonlinear equation numerically to minimize quantile_error

  return(sol)

}

# Loose, data-driven box for the BGEV parameters (mu, sigma, xi, delta). It only
# needs to comfortably contain plausible values and seed the perturbed starts;
# loose bounds act like starting values, not tight constraints. Returns a list
# with `lower` and `upper`.
bgev_start_box <- function(x, k = 3) {
  center <- stats::median(x)
  spread <- stats::mad(x)
  if (!is.finite(spread) || spread <= 0) spread <- stats::sd(x)
  if (!is.finite(spread) || spread <= 0) spread <- 1

  # mu around the median; sigma kept positive on a multiplicative scale;
  # xi over a loose fixed range; delta over a positive range, since estimation
  # restricts delta > 0 (see bgev_mle).
  list(lower = c(center - k * spread, spread / 10, -1, 0.05),
       upper = c(center + k * spread, spread * 10,  1, 5))
}


# Draw `n` uniform starting points inside a box from bgev_start_box().
bgev_sample_starts <- function(n, box) {
  if (n <= 0) return(matrix(numeric(0), ncol = length(box$lower)))
  p <- length(box$lower)
  out <- matrix(NA_real_, nrow = n, ncol = p)
  for (j in seq_len(p))
    out[, j] <- stats::runif(n, box$lower[j], box$upper[j])
  out
}


#' Profile log-likelihood for a BGEV parameter
#'
#' Computes the profile log-likelihood of one BGEV parameter over a grid,
#' maximising over the remaining three at each grid point (Nelder-Mead). This is
#' the diagnostic used in Otiniano et al. (2023) to check whether a fitted
#' optimum is global and the parameter is well identified.
#'
#' @param x Numeric vector of observations.
#' @param par Vector \code{c(mu, sigma, xi, delta)}, e.g. \code{bgev_mle(x)$par}.
#' @param which Index (1-4) or name of the parameter to profile.
#' @param span Half-width of the grid, as a fraction of the parameter value.
#' @param n Number of grid points.
#' @param plot Logical; if \code{TRUE}, plot the profile curve.
#'
#' @return A data frame with the grid values of the profiled parameter and the
#'   profile log-likelihood.
#'
#' @author Thiago do Rego Sousa
#'
#' @export
bgev_profile_likelihood <- function(x, par, which, span = 0.5, n = 41, plot = TRUE) {
  names_par <- c("mu", "sigma", "xi", "delta")
  j <- if (is.character(which)) match(which, names_par) else which
  center <- as.numeric(par[j])
  width <- max(abs(center) * span, span)
  grid <- seq(center - width, center + width, length.out = n)
  others <- setdiff(seq_len(4), j)

  prof <- vapply(grid, function(val) {
    obj <- function(rest) {
      pars <- numeric(4)
      pars[j] <- val
      pars[others] <- rest
      bgev_negative_log_likelihood(x, pars)
    }
    -stats::optim(par = as.numeric(par[others]), fn = obj, method = "Nelder-Mead")$value
  }, numeric(1))

  out <- data.frame(value = grid, loglik = prof)
  names(out)[1] <- names_par[j]
  if (plot) {
    plot(out[[1]], out$loglik, type = "l", xlab = names_par[j],
         ylab = "profile log-likelihood")
    graphics::abline(v = center, lty = 2)
  }
  out
}







