test_that("distCheck consistency holds for bgev with xi != 0 and delta != 0", {
  mu <- 0; sigma <- 1; xi <- -0.3; delta <- 1
  support <- bgev_support(mu, sigma, xi, delta)

  set.seed(42)
  res <- distCheck(fun = "bgev", n = 2000,
                    support.lower = support[1], support.upper = support[2],
                    subdivisions = 5000, mu = mu, sigma = sigma, xi = xi, delta = delta,
                    print.result = FALSE)

  expect_true(res$test1.density$error.check)
  expect_true(res$test2.quantile.cdf$error.check)
  expect_true(res$test3.mean.var$error.check)
})

test_that("density stays finite just inside the upper support boundary", {
  mu <- 0; sigma <- 1; xi <- -0.3; delta <- 1
  support <- bgev_support(mu, sigma, xi, delta)

  boundary_density <- dbgev(support[2] - 1e-10, mu = mu, sigma = sigma, xi = xi, delta = delta)

  expect_true(is.finite(boundary_density))
})
