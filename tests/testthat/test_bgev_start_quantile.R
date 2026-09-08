test_that("bgev_start_using_quantiles returns vector of size 4 for arbitrary input", {
  set.seed(1)
  x <- rnorm(200)
  res = bgev_start_using_quantiles(x)
  expect_equal(length(res), 4)
})
