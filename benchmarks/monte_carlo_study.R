source("R/bgev_domain.R")
source("R/bgev_distribution.R")
source("R/bgev_estimation.R")

library(SimDesign)

Design = SimDesign::createDesign(
  n = c(100, 500),
  mu = c(0),
  sigma = c(1),
  xi = c(-1,0,1),
  delta = c(0,1,9)
)

Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rbgev(n = n, mu = mu, sigma = sigma, xi = xi, delta = delta)  )
  dat
} 

Analyse <- function(condition, dat, fixed_objects) {

  ret_error <- rep(NA_real_, 4)
  names(ret_error) <- c("mu", "sigma", "xi", "delta")
  est <- tryCatch(bgev_mle(dat, control = DEoptim::DEoptim.control(itermax = 100, NP = 100, trace = FALSE)), error = function(e) NULL)
  if(is.null(est))
    return(ret_error)
  if (is.null(est$optim$bestmem) ||
      length(est$optim$bestmem) != 4 ||
      any(!is.finite(est$optim$bestmem))) {
    return(ret_error)
  }
  ret <- as.numeric(est$optim$bestmem)
  names(ret) <- c("mu", "sigma", "xi", "delta")
  return(ret)
}


Summarise <- function(condition, results, fixed_objects) {
  # assuming your Design object columns match these names
  true_mu <- condition$mu
  true_sigma <- condition$sigma
  true_xi <- condition$xi
  true_delta <- condition$delta
  
  # Return a named vector of the summary statistics (bias and RMSE)
  ret <- c(
    bias_mu = bias(results[, "mu"], parameter = true_mu),
    bias_sigma = bias(results[, "sigma"], parameter = true_sigma),
    bias_xi = bias(results[, "xi"], parameter = true_xi),
    bias_delta = bias(results[, "delta"], parameter = true_delta),
    RMSE_mu = RMSE(results[, "mu"], parameter = true_mu),
    RMSE_sigma = RMSE(results[, "sigma"], parameter = true_sigma),
    RMSE_xi = RMSE(results[, "xi"], parameter = true_xi),
    RMSE_delta = RMSE(results[, "delta"], parameter = true_delta) 
  )
  return(ret)
}

monte_carlo_results <- SimDesign::runSimulation(design=Design, replications=500,
                                  generate=Generate, analyse=Analyse, summarise=Summarise,
                                  progress = FALSE, verbose = FALSE, store_results = TRUE, 
                                  parallel = TRUE, ncores = 7, save_results = TRUE)

saveRDS(object = monte_carlo_results, file = "benchmarks/itermax_100_NP_100_replications_500.rds")
