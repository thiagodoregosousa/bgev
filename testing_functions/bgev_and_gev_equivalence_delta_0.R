#########
# TESTING EQUIVALENCE OF BGEV WITH DELTA = 0 AND GEV FROM ENVSTATS
##########
# new strategy based on the fact that when delta = 0, bgev equals gev with parametesr mu sigma and xi, 

# plot both densities to see if this is really true 
library(EnvStats)
x = seq(-5,5,0.01)
delta = 0.1
mu = 0
sigma = runif(1,0.5,2)
shape = runif(1,-2,2)
y_gev = dgevd(x, location = mu, scale = sigma, shape = -shape)
y_bgev = bgev::dbgev(x, mu = mu, sigma = sigma, xi = shape, delta = delta)
plot(x,y_gev, lwd = 10, lty = 2)
lines(x,y_bgev, lwd = 2, col = "blue")
y = dnorm(x,mean = pars_norm$estimate[1], sd = pars_norm$estimate[2])
hist(dt, breaks = 10, freq = FALSE)
