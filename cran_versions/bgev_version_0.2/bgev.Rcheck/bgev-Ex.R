pkgname <- "bgev"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bgev')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bgev-package")
### * bgev-package

flush(stderr()); flush(stdout())

### Name: bgev-package
### Title: Bimodal GEV Distribution with Location Parameter
### Aliases: bgev-package bgev dbgev pbgev qbgev rbgev
### Keywords: bgev gev bimodal distribution

### ** Examples

# generate 100 values distributed according to a bimodal GEV
x = rbgev(50, mu = 0.2, sigma = 1, xi = 0.5, delta = 0.2) 
# estimate the bimodal GEV parameters using the generated data
bgev.mle(x)



cleanEx()
nameEx("bgev.mle")
### * bgev.mle

flush(stderr()); flush(stdout())

### Name: bgev.mle
### Title: Parameter estimation of bimodal GEV distribution based on real
###   data that appears to be bimodal.
### Aliases: bgev.mle
### Keywords: bgev gev bimodal mle estimation

### ** Examples

# generate 100 values distributed according to a bimodal GEV
x = rbgev(50, mu = 0.2, sigma = 1, xi = 0.5, delta = 0.2) 
# estimate the bimodal GEV parameters using the generated data
bgev.mle(x)



cleanEx()
nameEx("bgev.support")
### * bgev.support

flush(stderr()); flush(stdout())

### Name: bgev.support
### Title: Support of the bimodal GEV distribution
### Aliases: bgev.support
### Keywords: bgev gev bimodal support

### ** Examples

# Computes the support of a specific bimodal GEV distribution
support = bgev.support(mu=1, sigma=10, xi=0.3, delta=2)



cleanEx()
nameEx("distCheck")
### * distCheck

flush(stderr()); flush(stdout())

### Name: distCheck
### Title: Check implementation of a distribution in R.
### Aliases: distCheck
### Keywords: distribution test implementation

### ** Examples

# generate random values for the parameters and test the 
# bimodal gev distribution implementation
set.seed(1)
mu = runif(1,-2,2)
sigma = runif(1,0.1,2)
xi = runif(1,0.3,0.9) * sign(runif(1,-1,1))
delta = 1#runif(1,-0.6,2)
support = bgev.support(mu, sigma, xi, delta)
var.exists = ( xi != 0) & (xi < (delta + 1)/2)
ret = distCheck(fun="bgev", n = 2000, 
      support.lower = support[1], support.upper = support[2], 
      subdivisions = 5000, mu = mu, sigma = sigma, xi = xi, 
      delta = delta, var.exists = var.exists, print.result = TRUE)



cleanEx()
nameEx("likbgev")
### * likbgev

flush(stderr()); flush(stdout())

### Name: likbgev
### Title: Log likelihood function for the bimodal GEV distribution.
### Aliases: likbgev
### Keywords: bgev gev bimodal likelihood

### ** Examples

# get random points from bimodal GEV
y = rbgev(100, mu = 1, sigma = 1, xi = 0.3, delta = 2)

# compute log-likelihood
likbgev (y, theta = c(1, 1, 0.3, 2)) 



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
