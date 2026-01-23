gev.support = function(mu, sigma, xi){
  xmax = 50
  support.lower = -xmax
  support.upper = xmax
  if( xi < 0)
    support.lower = mu + sigma/xi
  if( xi > 0)
    support.upper = mu + sigma/xi
  
  return(c(support.lower,support.upper))
}
gev.fExtremes.support = function(mu, sigma, xi){
  xmax = 50 
  support.lower = -xmax
  support.upper = xmax
  if( xi > 0)
    support.lower = mu - (sigma)/xi #(mu - sigma)/xi
  if( xi < 0)
    support.upper = mu - (sigma)/xi #(mu - sigma)/xi
  
  return(c(support.lower,support.upper))
}