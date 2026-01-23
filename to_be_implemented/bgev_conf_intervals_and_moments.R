# FUNCTIONS FOR FUTURE VERSIONS OF THE PACKAGE 
library(expint)



# auxiliar function
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )}


ebgevBOOTSTRAP_FUTURE <- function( x, brep, mu = 1, sigma = 1, xi = 0.3, delta = 2,conf.level = 0.95){
  # LOOK AT LIKELIHOOD, I CHANGED TO RETURN LOG LIKE INSTEAD OF - LOG LIKE. 
  
  # Description:
  #   Estimate the parameters of a BGEV using MLE and optionally construct a confidence interval for
  #   the parameters.
  #   Parameters: x a numeric vector of a BGEV distribution sample;
  #   brep > 0 number of bootstrap replicas;
  #   mu in R; sigma > 0; xi in R ;  delta > -1;
  #   conf.level is the confidence level for the confidence interval.
  
  # FUNCTION:
  # Error treatment of input parameters
  if(brep<0){
    stop("Number of breplicas must be greater than 0.")}
  if(sigma <= 0  || delta <= -1 ){
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")}
  
  # Create auxiliary variables:
  muhat      = rep(NA, times = brep)
  sigmahat   = rep(NA, times = brep)
  xihat      = rep(NA, times = brep)
  deltahat   = rep(NA, times = brep)
  
  starts=c(mu, sigma, xi, delta)
  n<-length(x)
  
  # Resampling with reposition
  for(k in (1:brep)){
    i <- sample(1:n, size = n, replace = TRUE)
    Z <- x[i]
    while(show_condition(suppressWarnings(optim(par= starts, fn = likbgev, y=Z, method="BFGS")))[1]=="error"){
      Z <- x[i]
    }
    
    esti <- optim(par= starts, fn = likbgev, y=Z, method="BFGS")
    
    
    muhat[k]       <- esti$par[1]
    sigmahat[k]    <- esti$par[2]
    xihat[k]       <- esti$par[3]
    deltahat[k]    <- esti$par[4]
  }
  
  estimate_mu     <-   mean(muhat,na.rm=TRUE)
  estimate_sigma  <-   mean(sigmahat,na.rm=TRUE)
  estimate_xi     <-   mean(xihat,na.rm = TRUE)
  estimate_delta  <-   mean(deltahat,na.rm=TRUE)
  
  
  CI_infmu<- estimate_mu-qnorm((1-conf.level/2),mean=0,sd=1)*sd(muhat,na.rm = TRUE)
  CI_supmu<-estimate_mu+qnorm((1-conf.level/2),mean=0,sd=1)*sd(muhat,na.rm = TRUE)
  
  CI_infsigma<- estimate_sigma-qnorm((1-conf.level/2),mean=0,sd=1)*sd(sigmahat,na.rm = TRUE)
  CI_supsigma<-estimate_sigma+qnorm((1-conf.level/2),mean=0,sd=1)*sd(sigmahat,na.rm = TRUE)
  
  CI_infxi<- estimate_xi-qnorm((1-conf.level/2),mean=0,sd=1)*sd(xihat,na.rm = TRUE)
  CI_supxi<-estimate_xi+qnorm((1-conf.level/2),mean=0,sd=1)*sd(xihat,na.rm = TRUE)
  
  CI_infdelta<- estimate_delta-qnorm((1-conf.level/2),mean=0,sd=1)*sd(deltahat,na.rm = TRUE)
  CI_supdelta<-estimate_delta+qnorm((1-conf.level/2),mean=0,sd=1)*sd(deltahat,na.rm = TRUE)
  
  cimu<-paste("(",CI_infmu,";",CI_supmu,")")
  cisigma<-paste("(",CI_infsigma,";",CI_supsigma,")")
  cixi<-paste("(",CI_infxi,";",CI_supxi,")")
  cidelta<-paste("(",CI_infdelta,";",CI_supdelta,")")
  
  tab <- as.table(rbind(c(estimate_mu , estimate_sigma , estimate_xi ,estimate_delta ), 
                        c(cimu, cisigma, cixi,cidelta)))
  dimnames(tab) <- list(Estimative = c("Ponctual", "Conf. Interval"),
                        parameter = c("mu", "sigma", "xi","delta"))
  print(tab)
  
}


smallgamma <- function(a,x){
  gammainc(a,0) - gammainc(a,x)
}



#----------------------------------------------------------------------
bgev.mean <- function(mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
  # Description:
  # Compute k-th moment E(X) for the  bimodal generalized extreme value distribution.
  # Reference: Cira EG Otiniano et al (2021). A Bimodal Model for Extremes Data.
  # Department of Statistics, University of Brasılia, Darcy Ribeiro,
  # Brasilia, 70910-900, DF, Brazil.
  # Department of Statistics, Federal University of Rio Grande do
  # Norte, Natal, 59078-970, RN, Brazil.
  # Parameters: mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
             sigma <= 0  || delta <= -1")
  
  
  if(delta == 0){
    return( (sigma/xi) * gamma(1-xi) + mu )
  }
  
  # Compute auxiliary variables
  mean.xi.positive  = (-1)^( (delta + 2) / (delta + 1) ) * (sigma/xi)^(1/(delta + 1)) *
    ( smallgamma(1-xi,1) - smallgamma(1,1) ) + 
    (sigma/xi)^(1/(delta + 1)) *  ( gammainc(1-xi,1) - gammainc(1,1) )
  # Return Value
  mean.xi.positive
}
#----------------------------------------------------------------------
