source("R/bgev_domain.R")
source("R/bgev_distribution.R")


log(x) = -1 
# quantile at y = exp(-1) = mu for bgev 

mu = 5 # runif(1,-1000,1000)
sigma =  runif(1,0.01,10)
xi = 0 # runif(1,-10,10)
delta = runif(1,-0.99,10)
qbgev(p = exp(-1), mu = mu, sigma = sigma, xi = xi, delta = delta) == mu


q1 = qbgev(p = exp(-exp(2)), mu = mu, sigma = sigma, xi = xi, delta = delta)
q2 = qbgev(p = exp(-exp(1)), mu = mu, sigma = sigma, xi = xi, delta = delta)
round(1/log(q1/q2, 2) - 1 - delta, 8) == 0
round(sigma == (-q2)^(delta+1), 8) == 0

mc_estimator <- function(n, R, mu, sigma, xi, delta) {
  sigma_constant = 0.1
  est <- matrix(NA, R, 3)
  colnames(est) <- c("mu_hat","sigma_hat","delta_hat")
  
  for(i in 1:R){
    
    x = rbgev(n = n, mu = mu, sigma = sigma, xi = xi, delta = delta)
    
    mu_hat = quantile(x, exp(-1))
    
    x_new = x - mu_hat
    
    q1 = quantile(x_new, exp(-exp(sigma_constant)))
    q2 = quantile(x_new, exp(-exp(1)))
    
    delta_hat = 1/log(q1/q2, sigma_constant) - 1
    sigma_hat = (-q2)^(delta_hat + 1)
    
    est[i,] = c(mu_hat, sigma_hat, delta_hat)
  }
  
  true_vals <- c(mu, sigma, delta)
  
  bias <- colMeans(est) - true_vals
  sd_est <- apply(est, 2, sd)
  
  list(
    estimates_mean = colMeans(est),
    bias = bias,
    sd = sd_est
  )
}

mc_estimator(100000, R = 100, mu = 3, sigma = 2, delta = 2, xi = 0)




xi = seq(-3,3,0.01)
plot(xi, 2^(-xi), type = 'l')
abline(h = 1)





library(nleqslv)
#' Starting values for BGEV distribution
#' 
#' @param x Numeric vector of observations.
#' @return starts_DEoptim start values for the estimation of (mu,sigma,xi,delta)
#' @author Thiago do Rego Sousa
bgev_start_using_quantiles = function(x){
  quantile_error <- function(theta){
    mu  <- theta[1]
    sigma <- theta[2]
    xi  <- theta[3]
    delta <- theta[4]
    
    p <- c(0.1,0.3,0.6,0.9)
    q_emp <- quantile(x, p)
    q_model <- qbgev(p, mu, sigma, xi, delta)
    
    return(q_model - q_emp)
  }
  
  start <- c(median(x), mad(x), 0.1, 1) # starting values for mu, sd, xi, delta
  sol <- nleqslv(start, quantile_error)$x # solve nonlinear equation numerically to minimize quantile_error
  
  return(sol)
  
}




bgev_start_region = function(x, n_samples = 200,     # LHS exploration points
                             n_params = 4,      # dimension
                             keep_frac =  0.15,    # keep best 15%
                             expand_factor =  0.10){  # expand refined box by 10%){
  delta = c(-0.99, 10)
  xi = c(-10,10)
  #q_10 = quantile(x, probs = 0.01)
  #q_90 = quantile(x, probs = 0.99)
  mu = c(min(x), max(x))  #c(q_10, q_90)
  sigma = c(0.1, (max(x) - min(x))*10)
  
  lower = c(mu[1],sigma[1],xi[1],delta[1])
  upper = c(mu[2],sigma[2],xi[2],delta[2])
  
  bgev_log_likelihood_negative = function(pars) {
    if(!bgev_valid_params(pars[1], pars[2], pars[3], pars[4]))
      return(1e+100)
    val = -bgev_log_likelihood(x, pars)
    if (!is.finite(val)) 
      return(1e+100)
    else return(val)
  }
  
  lhs_unit <- randomLHS(n_samples, n_params)
  
  samples <- matrix(NA, nrow = n_samples, ncol = n_params)
  
  for (i in 1:n_params) {
    samples[, i] <- lower[i] + lhs_unit[, i] * (upper[i] - lower[i])
  }
  
  values <- apply(samples, 1, bgev_log_likelihood_negative)
  
  
  ############################################
  # 5) SELECT BEST REGION
  ############################################
  
  n_keep <- ceiling(n_samples * keep_frac)
  best_idx <- order(values)[1:n_keep]
  best_points <- samples[best_idx, ]
  
  # Quantile-based refined region (robust)
  new_lower <- apply(best_points, 2, quantile, probs = 0.10)
  new_upper <- apply(best_points, 2, quantile, probs = 0.90)
  
  ############################################
  # 6) EXPAND REGION SLIGHTLY (SAFETY MARGIN)
  ############################################
  
  range_width <- new_upper - new_lower
  
  new_lower <- new_lower - expand_factor * range_width
  new_upper <- new_upper + expand_factor * range_width
  
  # Ensure region stays within original bounds
  new_lower <- pmax(new_lower, lower)
  new_upper <- pmin(new_upper, upper)
  
  # return 
  bounds = list()
  bounds$lower = lower
  bounds$upper = upper
  return(bounds)
  
  
  
}



total <- 100
count <- 0

for (i in 1:total) {
  
  # Generate true parameters
  mu    <- runif(1, -1000, 100)
  sigma <- runif(1, 0.01, 1000)
  xi    <- runif(1, -10, 10)
  delta <- runif(1, -0.99, 10)
  
  # Simulate data with true parameters
  x <- rbgev(n = 1000, mu = mu, sigma = sigma, xi = xi, delta = delta)
  
  # Get estimated starting region
  start_result <- bgev_start_region(x)
  
  lw <- start_result$lower
  up <- start_result$upper
  
  # Check if true parameters lie inside bounds
  inside = FALSE
  inside <- (
    lw[1] <= mu    && mu      <= up[1] &&
      lw[2] <= sigma && sigma <= up[2] &&
      lw[3] <= xi    && xi    <= up[3] &&
      lw[4] <= delta && delta <= up[4]
  )
  
  if(inside == TRUE)
    print(c(i,count))
  
  if (inside) {
    count <- count + 1
  }
}

# Percentage of times bounds contain true parameters
percentage <- 100 * count / total

cat("Coverage percentage:", percentage, "%\n")

