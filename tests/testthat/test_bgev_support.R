test_that("bgev_support matches the closed-form GEV support at delta = 0", {
  # Frechet-type (xi > 0): bounded below at mu - sigma/xi
  expect_equal(bgev_support(mu = 0, sigma = 1, xi = 0.5, delta = 0)[1], 0 - 1/0.5)

  # Weibull-type (xi < 0): bounded above at mu - sigma/xi
  expect_equal(bgev_support(mu = 0, sigma = 1, xi = -0.5, delta = 0)[2], 0 - 1/(-0.5))

  # Gumbel-type (xi = 0): unbounded on both sides
  expect_equal(bgev_support(mu = 0, sigma = 1, xi = 0, delta = 0), c(-Inf, Inf))
})
