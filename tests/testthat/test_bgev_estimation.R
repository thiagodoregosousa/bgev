test_that("asdf", {

  x = rbgev(500, mu = -5, sigma = 0.1, xi = 3, delta = 0)
  est = bgev_mle(x)

   #expect_equal(density_bgev, density_gev)
})