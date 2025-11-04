################################################################################
# FUNCTION:    Bimodal GEV (generalized extreme value) distribution 
#              proposed in Cira EG Otiniano, Bianca S Paiva, Roberto 
#              Vila and Marcelo Bourguignon (2021)
#  dbgev       Density for the bimodal generalized extreme value distribution.
#  qbgev       Quantile function for the bimodal GEV distribution. 
#  rbgev       Random generation for the bimodal GEV distribution.
#  likbgev     maximum likelihood (ML) estimators for the parameters of a BGEV
#              distribution
#  ebgev       Estimate the parameters of a BGEV and optionally construct 
#              a confidence interval for the parameters.
################################################################################


#library(EnvStats)
#library(DEoptim)



#-------------------------------------------------------------------------------
dbgev <- function(y, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
  # DESCRIPTION:
  # Compute the density for the bimodal generalized extreme value distribution.
  # Parameters: y in R; mu in R; sigma > 0; xi in R; delta > -1;
    
  # FUNCTION:
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  
  # Compute auxiliary variables:
  T      <- (y-mu)*(abs(y-mu)^delta)
  derivate_T <- (delta + 1)*(abs(y-mu)^delta)
  
  # Compute density points
  pdf    <- dgevd(T, 0, scale=sigma, shape=-xi)*derivate_T # changed shape to -xi due to reparametrization
  
  # Return Value
  pdf
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
pbgev <- function(y, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
  # DESCRIPTION:
  # Compute the distribution function for the bimodal generalized extreme value distribution.
  # Parameters: y in R; mu in R; sigma > 0; xi in R; delta > -1;
  
  # FUNCTION:
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  
  # Compute auxiliary variables:
  Ti      <- (y-mu)*(abs(y-mu)^delta)
  # Compute 
  cdf    <- pgevd(Ti, location=0, scale=sigma, shape=-xi)   # changed shape to -xi due to reparametrization
  # Return Value
  return(cdf)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
qbgev   <- function(p, mu = 1, sigma = 1, xi = 0.3, delta = 2){
  # DESCRIPTION:
  # Compute the quantile for the Bimodal GEV distribution
  # Parameters:  p in [0;1];  mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  
  # Compute distribution points according to their sign
  quantile <- sign(qgevd(p, 0, sigma, -xi))*(abs(qgevd(p, 0, sigma, -xi)))^(1/(delta + 1)) + mu    # changed shape to -xi due to reparametrization
  # Return Value
  return(quantile)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
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
  U <- runif(n)
  # Compute random numbers
  rnumber <- qbgev(U, mu, sigma, xi, delta)
  # Return Value
  return(rnumber)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
likbgev <- function(y, theta = c(1, 1, 0.3, 2)){
  # DESCRIPTION:
  # log-likelihood function for the parameters of a BGEV distribution
  # Parameters: y in R; theta: vector with mu, sigma, xi and delta, respectively.
  # mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  mu      <- theta[1]
  sigma   <- theta[2]
  xi      <- theta[3]
  delta   <- theta[4]
  
  # Error treatment of input parameters
  if(length(theta)!=4){
    stop("vector of parameters needs to be of length 4.")}
  if(sigma <= 0  || delta <= -1 ){
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")}
  
  # Compute auxiliary variables:
  T      <- (y-mu)*(abs(y-mu)^delta)
  derivate_T <- (delta + 1)*(abs(y-mu)^delta)
  # Compute density points
  dbgevd <- dgevd(x = T, location = mu, scale = sigma, shape = -xi)*derivate_T
  # Log:
  logl <- sum(log(dbgev(y,mu, sigma, xi, delta)))
  return(logl)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
bgev.mle <- function(x, lower = c(-3,0.1,-3,-0.9), upper = c(3,3,3,3), 
                     itermax = 50){
  # DESCRIPTION:
  # log-likelihood function for the parameters of a BGEV distribution
  # Parameters: x: data for estimating the parameters of the bimodal GEV
  # lower and uppper: lower and upper bounds to search for the true parameter
  # itermax: maximum number of interations when finding good starting 
  # values using DEOptim
  
  # FUNCTION:
  # get reasonable starting values using a genetic algorithm 
  lower = c(-3,0.1,-3,-0.9)
  upper = c(3,3,3,3)
  likbgev2 = function(theta){
    # function to be minimized
    val = -likbgev(x, theta)
    if (val == Inf)
      return(1e100)
    else
      return(val)
  }
  # theta = mu,sigma,xi,delta, # DEoptim minimizes
  starts = DEoptim(fn = likbgev2, lower, upper, control = DEoptim.control(itermax = itermax, trace = FALSE)) 

  # second optimization step starting at the value returned by deoptim. Recall that optim minimizes
  esti <- optim(par = starts$optim$bestmem, fn = likbgev2, method="L-BFGS-B", lower = lower, upper = upper)
  
  # return
  esti
    
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
bgev.support = function(mu, sigma, xi, delta){
  # DESCRIPTION:
  # Computes the support in R of the specified bimodal GEV distribution.
  # In the functions that follow, 'maxVal' should be ideally +Inf, 
  # but you can use a value to avoid numerical integration errors when using 
  # functions that depend on the density or distribution function. 
  
  # FUNCTION:  
  maxVal = Inf
  support.lower = -maxVal
  support.upper = maxVal
  if( xi > 0)
    support.lower = mu - (sigma/xi)^(1/(delta+1))
  if( xi < 0)
    support.upper = mu + abs(sigma/xi)^(1/(delta+1))
  
  # return
  return(c(support.lower,support.upper))
}
#-------------------------------------------------------------------------------



