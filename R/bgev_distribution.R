#' Bimodal GEV (generalized extreme value) distribution
#'
#' Functions to compute the density, distribution function, quantile function,
#' and to generate random variates for the bgebv (bimodal generalized extreme value)
#'
#' This distribution corresponds was proposed by in Cira EG Otiniano, 
#' Bianca S Paiva, Roberto Vila and Marcelo Bourguignon (2021)
#'
#' @name bgev
#' @rdname bgev
#' @aliases bgev dbgev pbgev qbgev rbgev bgev_support
#'
#' @param x Numeric vector of values for calculating density. 
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations for random generation.
#' @param mu location parameter
#' @param sigma scale parameter (sigma > 0)
#' @param xi shape parameter in R
#' @param delta shape parameter (delta > -1)
#' @param log Logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return
#' \item{dbgev}{density values}
#' \item{pbgev}{distribution function values}
#' \item{qbgev}{quantile function values}
#' \item{rbgev}{random variates}
#' 
#' @note
#' BGEV distribution is equivalent to the GEV distribution when \code{delta = 0}. 
#' When comparing BGEV with GEV from package EnvStats, the \code{shape} parameter of GEV
#' is changed to \code{-xi} due to reparametrization
#' 
#' 
#'
#' @references
#' Otiniano, Cira E. G., et al. (2023).
#' \emph{A bimodal model for extremes data}.
#' Environmental and Ecological Statistics, 1–28.
#' \doi{10.1007/s10651-023-00566-7}
#'
#' @author Thiago do Rego Sousa
#'
#' @examples
#' par(mfrow = c(2, 2))
#' set.seed(1000)
#' r <- rbgev(n = 1000)
#' plot(r, type = "l", main = "BGEV Random Values")
#'
#' hist(r, probability = TRUE, border = "white", ylim = c(0,1))
#' x <- seq(min(r), max(r), length = 201)
#' lines(x, dbgev(x), lwd = 2)
#'
#' plot(sort(r), (1:1000)/1000, main = "Probability", ylab = "Probability")
#' lines(x, pbgev(x), lwd = 2)
#'
#' round(qbgev(pbgev(q = seq(0, 3, by = 0.1)), 6),2)
#' @export
dbgev <- function(x, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
  
  if(!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma must be > 0 and delta must be > -1")

  support = bgev_support(mu, sigma, xi, delta)
  support_lower = support[1]
  support_upper = support[2]

  if( xi > 0)
    out_of_support_indexes = which(x < support_lower)
  if( xi < 0)
    out_of_support_indexes = which(x > support_upper)

  pdf = rep(0, length(x))
  inside_support_indexes = setdiff(1:length(x), out_of_support_indexes)

  x_in = x[inside_support_indexes]

  # Compute auxiliary variables:
  T      <- (x_in-mu)*(abs(x_in-mu)^delta)
  derivate_T <- (delta + 1)*(abs(x_in-mu)^delta)

  # Compute density points
  pdf_x_in    <- EnvStats::dgevd(T, 0, scale=sigma, shape=-xi)*derivate_T # changed shape to -xi due to reparametrization
  
  # Return Value
  pdf[inside_support_indexes] = pdf_x_in
  pdf
}


#' @export
pbgev <- function(q, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 

   if(!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma must be > 0 and delta must be > -1")
  
  # Compute auxiliary variables:
  Ti      <- (q-mu)*(abs(q-mu)^delta)
  cdf    <- EnvStats::pgevd(Ti, location=0, scale=sigma, shape=-xi)   # changed shape to -xi due to reparametrization
  return(cdf)
}


#' @export
qbgev   <- function(p, mu = 1, sigma = 1, xi = 0.3, delta = 2){
  if(!bgev_valid_params(mu, sigma, xi, delta))
    stop("Invalid parameters: sigma must be > 0 and delta must be > -1")
  
  # Compute distribution points according to their sign
  quantile <- sign(EnvStats::qgevd(p, 0, sigma, -xi))*(abs(EnvStats::qgevd(p, 0, sigma, -xi)))^(1/(delta + 1)) + mu    # changed shape to -xi due to reparametrization
  return(quantile)
}


#' @export 
rbgev <- function(n, mu = 1, sigma = 1, xi = 0.3, delta = 2){
  # DESCRIPTION:
  # random generator for the Bimodal GEV distribution.
  # Parameters:  n in {1,2,3,...}; mu in R; sigma > 0; xi in R ; delta > -1;
  
  # FUNCTION:
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  # Compute auxiliary variables:
  U <- stats::runif(n)
  # Compute random numbers
  rnumber <- qbgev(U, mu, sigma, xi, delta)
  # Return Value
  return(rnumber)
}