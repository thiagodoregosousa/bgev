#' Log-likelihood function for the BGEV distribution
#'
#' @param x Numeric vector of observations.
#' @param theta Vector of parameters (mu, sigma, xi, delta). See \link{bgev}.
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#'
#' @return The log-likelihood value.
#' @export 
bgev_log_likelihood <- function(x, theta = c(1, 1, 0.3, 2)){

  mu      <- theta[1]
  sigma   <- theta[2]
  xi      <- theta[3]
  delta   <- theta[4]
  
  if(!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma must be > 0 and delta must be > -1")
  
  # Log-likelihood:
  logl <- sum(log(dbgev(x,mu, sigma, xi, delta)))
  return(logl)
}

#' Maximum Likelihood Estimation for the BGEV distribution
#' 
#' @param x Numeric vector of observations.
#' @param deoptim.itermax Maximum number of iterations for the DEoptim algorithm to get for all parameters 
#' @param optim.method Optimization method to be used in \code{optim} for local optimization. See \code{method} argument in \code{?optim}.
#' @param start Optional vector of starting values for the parameters (mu, sigma, xi, delta).
#' @param lower Optional vector of lower bounds for the parameters (mu, sigma, xi, delta).
#' @param upper Optional vector of upper bounds for the parameters (mu, sigma, xi, delta).
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#' 
#' @note lower and upper should be provided together.
bgev_mle <- function (x, method_envstats = "mle", deoptim.itermax = 200, 
          optim.method = "L-BFGS-B", start = NULL, lower = NULL, upper = NULL, verbose = FALSE) 
{
  likbgev2 = function(theta) {
    val = -likbgev(x, theta)
    if (!is.finite(val)) 
      return(1e+100)
    else return(val)
  }

 # get lower and upper bounds if not provided 
 if (is.null(lower) | is.null(upper)) {
    fit_gev <- try(EnvStats::egevd(x = x, method = method_envstats))
    if (inherits(fit_gev, "try-error")) {
      if(verbose)
        message("I tried to get good starting values using egevd and it failed. Try using a different method for method_envstats, e.g, pwme. See the help of egevd for the available options.")
      fit_gev <- try(EnvStats::egevd(x = x, method = "pwme"))
      #return(NULL)
    }
    mu_start = fit_gev$parameters[1]
    sd_start = fit_gev$parameters[2]
    xi_start = -fit_gev$parameters[3]
    xi_min = xi_start/5
    xi_max = xi_start * 5
    if (xi_start < 0) {
      xi_max = xi_start/5
      xi_min = xi_start * 5
    }
    lower = c(mu_start - 2 * sd_start, sd_start/5, xi_min, 
              -0.99)
    upper = c(mu_start + 2 * sd_start, 5 * sd_start, xi_max, 
              20)
  }
  starts.DEoptim = DEoptim::DEoptim(fn = likbgev2, lower, upper, 
                                    control = DEoptim::DEoptim.control(itermax = deoptim.itermax, 
                                                                       trace = FALSE))
  esti <- stats::optim(par = starts.DEoptim$optim$bestmem, 
                       fn = likbgev2, method = optim.method, lower = lower, 
                       upper = upper)
  esti
}
