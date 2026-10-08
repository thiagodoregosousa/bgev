
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
