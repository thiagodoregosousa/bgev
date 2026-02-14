#' Log-likelihood function for the BGEV distribution
#'
#' @param x Numeric vector of observations.
#' @param pars Vector of parameters (mu, sigma, xi, delta). See \link{bgev}.
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#'
#' @return The log-likelihood value.
#' @export 
bgev_log_likelihood <- function(x, pars = c(1, 1, 0.3, 2)){

  mu      <- pars[1]
  sigma   <- pars[2]
  xi      <- pars[3]
  delta   <- pars[4]
  
  if(!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma must be > 0 and delta must be > -1")
  
  # Log-likelihood:
  logl <- sum(log(dbgev(x,mu, sigma, xi, delta)))
  return(logl)
}


#' Maximum Likelihood Estimation for the BGEV distribution
#' 
#' @param x Numeric vector of observations.
#' @param control List of type DEoptim::DEoptim.control (PUT LINK HERE)
#' @param lower Optional vector of lower bounds for the parameters (mu, sigma, xi, delta).
#' @param upper Optional vector of upper bounds for the parameters (mu, sigma, xi, delta).
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#' 
bgev_mle <- function (x, lower = c(-12,0.01,-0.99,-12), upper = c(12,12,12,12), 
                       control = DEoptim::DEoptim.control(itermax = 20, NP = 100, trace = FALSE), 
                       DEoptim_replicates = 5) 
{
  if (is.null(x) || anyNA(x)) {
    stop("`x` must be a non-null numeric vector with no missing values.", call. = FALSE)
  }
  
  bgev_log_likelihood_negative = function(pars) {
    if(!bgev_valid_params(pars[1], pars[2], pars[3], pars[4]))
      return(1e+100)
    val = -bgev_log_likelihood(x, pars)
    if (!is.finite(val)) 
      return(1e+100)
    else return(val)
  }

  fits <- replicate(DEoptim_replicates, {
    DEoptim::DEoptim(fn = bgev_log_likelihood_negative, control = control, lower = lower, upper = upper)
  }, simplify = FALSE)
  
  best <- fits[[ which.max(sapply(fits, function(f) -f$optim$bestval)) ]]
  
  best
}







