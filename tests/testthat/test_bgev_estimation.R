

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


test_that("bgev_mle validates x", {
  expect_error(bgev_mle(NULL))
  expect_error(bgev_mle("abc"))
  expect_error(bgev_mle(c(NA, 1, 2)))
})

test_that("bgev_mle depends on x", {
  set.seed(1)
  r1 <- bgev_mle(rnorm(100))
  set.seed(1)
  r2 <- bgev_mle(rnorm(100, mean = 3))
  expect_false(isTRUE(all.equal(r1$optim$bestmem, r2$optim$bestmem)))
})

test_that("estimates respect bounds", {
  set.seed(1)
  fit <- bgev_mle(rnorm(200))
  expect_true(all(fit$optim$bestmem >= fit$member$lower))
  expect_true(all(fit$optim$bestmem <= fit$member$upper))
})


test_that("results are reproducible with fixed seed", {
  set.seed(42)
  r1 <- bgev_mle(rnorm(200))
  
  set.seed(42)
  r2 <- bgev_mle(rnorm(200))
  
  expect_equal(r1$par, r2$par, tolerance = 1e-6)
})


test_that("optimization improves likelihood", {
  set.seed(1)
  x <- rnorm(200)
  
  fit <- bgev_mle(x)
  
  random_par <- runif(4, fit$member$lower, fit$member$upper)
  
  expect_gte(
    fit$optim$bestval,
    bgev_log_likelihood(x, random_par)
  )
})


test_that("DEoptim control affects result", {
  set.seed(1)
  x <- rnorm(200)
  r1 <- bgev_mle(x,
                 control = DEoptim::DEoptim.control(itermax = 5, NP = 100, trace = FALSE)
  )
  
  set.seed(1)
  x <- rnorm(200)
  r2 <- bgev_mle(x,
                 control = DEoptim::DEoptim.control(itermax = 10, NP = 100, trace = FALSE)
  )
  
  expect_false(isTRUE(all.equal(r1$optim$bestmem, r2$optim$bestmem)))
})



test_that("replicates improve or match solution", {
  set.seed(1)
  x <- rnorm(200)
  
  r1 <- bgev_mle(x, DEoptim_replicates = 1)
  r5 <- bgev_mle(x, DEoptim_replicates = 5)
  
  expect_lte(r5$optim$bestval, r1$optim$bestval)
})







