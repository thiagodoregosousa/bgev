

#' Compute the support of the BGEV distribution
#' 
#' Returns the lower and upper limits of the support of the BGEV distribution
#' 
#' @param mu location parameter
#' @param sigma scale parameter (sigma > 0)
#' @param xi shape parameter in R
#' @param delta shape parameter (delta > -1)
#' 
#' #' @author Thiago do Rego Sousa
#' 
#' @return A vector of length 2 with the lower and upper limits of the support
#' 
#' @details It returns values with \code{-Inf} or \code{Inf} when the support is unbounded. 
#' When the shape parameter \code{xi} is different from zero, the support 
#' is truncated either at the left or at the right side of the real. 
#' Considering the support is particularly useful to estimating momoments and 
#' to compute the likelihood function.
bgev_support = function(mu = 1, sigma = 1, xi = 0.3, delta = 2){
  support_lower = -Inf
  support_upper = Inf
  if( xi > 0)
    support_lower = mu - (sigma/xi)^(1/(delta+1))
  if( xi < 0)
    support_upper = mu + abs(sigma/xi)^(1/(delta+1))  
  return(c(support_lower,support_upper))
}


#' Validate BGEV parameters
#' 
#' Check if the provided parameters for the BGEV distribution are valid.
#' 
#' @param \link{bgev} Parameters of the BGEV distribution as in \link{bgev}
#' 
#' @author Thiago do Rego Sousa
#' 
#' @export 
bgev_valid_params = function(mu, sigma, xi, delta){
  if(sigma <= 0  || delta <= -1 )
    return(FALSE)
  return(TRUE)
}
