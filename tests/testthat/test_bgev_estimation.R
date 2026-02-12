



test_that("likelihood behaviour at invalid parameters", {
  # invalid sigma
  expect_error(bgev_log_likelihood(rbgev(100), pars = c(0,0,1,1)))
  
  # invalid delta
  expect_error(bgev_log_likelihood(rbgev(100), pars = c(0,1,1,-1.1)))
  
})


test_that("likelihood finite at truth", {

  x <- rbgev(n = 100, 1, 1, 1, 1)

  log_likelihood_x = bgev_log_likelihood(x, pars = c(1,1,1,1))

  expect_true(is.finite(log_likelihood_x))

  expect_true(log_likelihood_x != 1e99)

})


test_that("likelihood worsen when parameters are perturbed", {

  
  pars_true <- c(1, 1, 1, 1)
  x <- rbgev(n = 100, pars_true[1], pars_true[2], pars_true[3], pars_true[4])
  log_likelihood_pars_true = bgev_log_likelihood(x, pars = pars_true)

  # pertub each of the 4 parameters by 2
  for( i in 1:4){
    pars_true_perturbed = pars_true
    pars_true_perturbed[i] = pars_true_perturbed[i] + 2
    log_likelihood_pars_true_perturbed = bgev_log_likelihood(x, pars_true_perturbed)

    expect_gt(log_likelihood_pars_true, log_likelihood_pars_true_perturbed )
  }

})