#-------------------------------------------
# gev
# mu in R; sigma > 0; xi in R

source("BGEV_functions.R")
source("distCheck.R")
library(fBasics)

mu = runif(1,-2,2)
sigma = runif(1,0.1,2)
xi = runif(1,0.1,2) * sign(runif(1,-1,1))
#mu = 0
#sigma = 1
#xi = 0.5
support = gev.support(mu, sigma, xi)
#support = c(-50,50)
distCheck(fun="gevd", n = 2000, 
          support.lower = support[1], support.upper = support[2], subdivisions = 5000,
          location = mu, scale = sigma, shape = xi)
xval = seq(support[1],support[2],0.01)
xval = seq(-5,5,0.01)
yval = dgevd(xval, mu, sigma, xi)
plot(xval,yval,type = 'b')

# so I tried with the GEV as well and got the error as well. 
# So somehow the support values do not work well. 
fBasics::distCheck("gevd", location = mu, scale = sigma, shape = xi)





# gev test with fExtremes
mu = 0
sigma = 1
xi = 0.5
mu = runif(1,-2,2)
sigma = runif(1,0.1,2)
xi = runif(1,0.1,2) * sign(runif(1,-1,1))
support = gev.fExtremes.support(mu, sigma, xi)
distCheck(fun="gev", n = 2000, 
          support.lower = support[1], support.upper = support[2], subdivisions = 1000,
          mu = mu, beta = sigma, xi = xi)
c(xi,mu)
xval = seq(support[1],support[2],0.01)
xval = seq(-5,5,0.01)
yval = fExtremes::dgev(xval,  mu = mu, beta = sigma, xi = xi)
plot(xval,yval,type = 'b', main = round(support,2))

