test_that("likelihood behaviour at invalid parameters", {
  # invalid sigma
  expect_error(bgev_log_likelihood(rbgev(100), pars = c(0, 0, 1, 1)))

  # invalid delta
  expect_error(bgev_log_likelihood(rbgev(100), pars = c(0, 1, 1, -1.1)))
})


test_that("likelihood finite at truth", {
  x <- rbgev(n = 100, 1, 1, 1, 1)
  log_likelihood_x <- bgev_log_likelihood(x, pars = c(1, 1, 1, 1))

  expect_true(is.finite(log_likelihood_x))
})


test_that("likelihood worsen when parameters are perturbed", {
  pars_true <- c(1, 1, 1, 1)
  x <- rbgev(n = 100, pars_true[1], pars_true[2], pars_true[3], pars_true[4])
  log_likelihood_pars_true <- bgev_log_likelihood(x, pars = pars_true)

  # perturb each of the 4 parameters by 2
  for (i in 1:4) {
    pars_true_perturbed <- pars_true
    pars_true_perturbed[i] <- pars_true_perturbed[i] + 2
    log_likelihood_pars_true_perturbed <- bgev_log_likelihood(x, pars_true_perturbed)

    expect_gt(log_likelihood_pars_true, log_likelihood_pars_true_perturbed)
  }
})


test_that("bgev_mle validates x", {
  expect_error(bgev_mle(NULL))
  expect_error(bgev_mle("abc"))
  expect_error(bgev_mle(c(NA, 1, 2)))
})


test_that("bgev_mle returns a 4-parameter estimate with diagnostics", {
  set.seed(1)
  fit <- bgev_mle(rnorm(200))

  expect_equal(length(fit$par), 4)
  expect_equal(names(fit$par), c("mu", "sigma", "xi", "delta"))
  expect_true(is.finite(fit$loglik))
  expect_true(is.logical(fit$agree))
  expect_true(is.list(fit$boundary))
  expect_true(all(is.finite(fit$par)))
  # reparametrisation keeps the estimate in the valid domain
  expect_gt(fit$par[["sigma"]], 0)
  expect_gt(fit$par[["delta"]], -1)
})


test_that("bgev_mle reports optimum diagnostics and admissibility", {
  set.seed(1)
  # well-identified BGEV sample: gradient ~ 0 and a positive-definite Hessian
  x <- rbgev(n = 300, mu = 0, sigma = 1, xi = 0.2, delta = 1)
  fit <- bgev_mle(x)

  expect_true(is.list(fit$optimum))
  expect_true(all(c("grad_norm", "hessian_pd", "eig_ratio") %in% names(fit$optimum)))
  expect_lt(fit$optimum$grad_norm, 1)      # near-zero on the negative-loglik scale
  expect_true(fit$optimum$hessian_pd)      # genuine, well-identified maximum
  expect_true(fit$admissible)              # accepted via the PD-Hessian gate
})


test_that("bgev_mle restricts estimation to delta > 0", {
  set.seed(1)
  fit <- bgev_mle(rbgev(n = 300, mu = 0, sigma = 1, xi = 0.2, delta = 1))
  expect_gt(fit$par[["delta"]], 0)
})


test_that("bgev_mle warns on tied data only under the continuous density", {
  set.seed(1)
  x <- sample(10:40, 80, replace = TRUE)  # integer ties
  expect_warning(bgev_mle(x), "tied")
  # the grouped likelihood is the right model for discrete data: no tie warning
  expect_warning(bgev_mle(x, likelihood = "grouped_likelihood"), NA)
})


test_that("bgev_mle grouped likelihood fits integer data", {
  set.seed(1)
  x <- sample(10:40, 120, replace = TRUE)
  fit <- bgev_mle(x, likelihood = "grouped_likelihood", h = 1)

  expect_equal(fit$likelihood, "grouped_likelihood")
  expect_equal(length(fit$par), 4)
  expect_gt(fit$par[["delta"]], 0)
  expect_true(is.finite(fit$loglik))
})


test_that("bgev_mle rejects an unknown likelihood", {
  expect_error(bgev_mle(rnorm(50), likelihood = "nope"))
})


test_that("bgev_mle depends on x", {
  set.seed(1)
  r1 <- bgev_mle(rnorm(100))
  set.seed(1)
  r2 <- bgev_mle(rnorm(100, mean = 3))
  expect_false(isTRUE(all.equal(r1$par, r2$par)))
})


test_that("results are reproducible with fixed seed", {
  set.seed(42)
  r1 <- bgev_mle(rnorm(200))

  set.seed(42)
  r2 <- bgev_mle(rnorm(200))

  expect_equal(r1$par, r2$par, tolerance = 1e-6)
})


test_that("reported loglik matches the log-likelihood at the estimate", {
  set.seed(1)
  x <- rnorm(200)
  fit <- bgev_mle(x)

  expect_equal(fit$loglik, bgev_log_likelihood(x, fit$par), tolerance = 1e-6)
})


test_that("more starts reach at least as good an admissible optimum", {
  # on well-identified BGEV data both runs find the same interior maximum;
  # the acceptance gate makes loglik non-monotone on misspecified data, so this
  # invariant only holds where an admissible optimum exists.
  set.seed(1)
  x <- rbgev(n = 300, mu = 0, sigma = 1, xi = 0.2, delta = 1)

  set.seed(1)
  r1 <- bgev_mle(x, n_starts = 1)
  set.seed(1)
  r10 <- bgev_mle(x, n_starts = 10)

  expect_gte(r10$loglik, r1$loglik - 1e-2)
})


test_that("bgev_start_using_quantiles returns vector of size 4 for arbitrary input", {
  set.seed(1)
  x <- rnorm(200)
  res <- bgev_start_using_quantiles(x)

  expect_equal(length(res), 4)
})


test_that("bgev_profile_likelihood returns a grid of the requested size", {
  set.seed(1)
  x <- rnorm(200)
  fit <- bgev_mle(x)
  prof <- bgev_profile_likelihood(x, fit$par, which = "xi", n = 11, plot = FALSE)

  expect_equal(nrow(prof), 11)
  expect_equal(names(prof), c("xi", "loglik"))
})
