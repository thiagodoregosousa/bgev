source("R/bgev_functions.R")
library(SimDesign)

Design = SimDesign::createDesign(
  n = c(100, 500),
  mu = c(-5, 5),
  sigma = c(0.1, 5),
  xi = c(0.1, 3),
  delta = c(0, 3)
)

Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rbgev(n = n, mu = mu, sigma = sigma, xi = xi, delta = delta)  )
  dat
} 

Analyse <- function(condition, dat, fixed_objects) {
  est = bgev_mle(dat)
  if(est$value == 0)
    return(rep(NA,4))
  ret = as.vector(est$par)
  names(ret) = c("mu", "sigma", "xi", "delta")
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
Design = Design[7:10, ]
monte_carlo_results <- SimDesign::runSimulation(design=Design, replications=5,
                                  generate=Generate, analyse=Analyse, summarise=Summarise,
                                  progress = FALSE, verbose = FALSE)
t(monte_carlo_results)
