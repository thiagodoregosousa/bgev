source("R/bgev_domain.R")
source("R/bgev_distribution.R")
source("R/bgev_estimation.R")

library(SimDesign)

Design = SimDesign::createDesign(
  n = c(100, 1000),
  mu = c(-50,0,50),
  sigma = c(0.1,5,50),
  xi = c(-10,0,1,10),
  delta = c(-0.5,1,10)
)

Design = SimDesign::createDesign(
  n = c(1000),
  mu = c(20),
  sigma = c(20),
  xi = c(2),
  delta = c(2)
)


Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rbgev(n = n, mu = mu, sigma = sigma, xi = xi, delta = delta)  )
  dat
} 

Analyse <- function(condition, dat, fixed_objects) {
  
  ret_error <- rep(NA_real_, 4)
  names(ret_error) <- c("mu", "sigma", "xi", "delta")
  est <- tryCatch(bgev_start_using_quantiles(dat), error = function(e) NULL)
  if(is.null(est))
    return(ret_error)
  if (length(est) != 4) {
    return(ret_error)
  }
  ret <- est
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

monte_carlo_results <- SimDesign::runSimulation(design=Design, replications=5,
                                                generate=Generate, analyse=Analyse, summarise=Summarise,
                                                progress = FALSE, verbose = FALSE, store_results = TRUE, 
                                                parallel = TRUE, ncores = 7, save_results = TRUE)
t(monte_carlo_results)


#saveRDS(object = monte_carlo_results, file = "benchmarks/itermax_100_NP_100_replications_500.rds")

