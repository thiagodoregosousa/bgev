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

#' Starting values for BGEV distribution
#'
#' @param x Numeric vector of observations.
#' @return starts_DEoptim start values for the estimation of (mu,sigma,xi,delta)
#' @author Thiago do Rego Sousa
#' @importFrom stats mad median quantile sd
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
  sol <- nleqslv::nleqslv(start, quantile_error)$x # solve nonlinear equation numerically to minimize quantile_error

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

  lhs_unit <- lhs::randomLHS(n_samples, n_params)

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
