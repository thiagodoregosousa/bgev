mu = runif(1,-2,2)
sigma = runif(1,0.1,2)
xi = runif(1,0.1,2) * sign(runif(1,-1,1))
delta = runif(1,-0.9,2)
support = bgev.support(mu, sigma, xi, delta)
x = seq(-10,10,0.01)
var.exists = ( xi != 0) & (xi < (delta + 1)/2)
y = dbgev(x, mu = mu, sigma = sigma, xi = xi, delta = delta)
plot(x,y, main = round(support,2))
abline(v = support)


x = seq(-1.1,-1.08,0.0001)
var.exists = ( xi != 0) & (xi < (delta + 1)/2)
y = dbgev(x, mu = mu, sigma = sigma, xi = xi, delta = delta)
plot(x,y)
abline(v = support)


dbgev(support[2]-1e-10, mu = mu, sigma = sigma, xi = xi, delta = delta)


distCheck(fun="bgev", n = 2000, 
          support.lower = -Inf, support.upper = support[2], subdivisions = 5000,
          mu = mu, sigma = sigma, xi = xi, delta = delta, var.exists = var.exists, robust = TRUE)
c(mu,sigma,xi,delta)
mu
