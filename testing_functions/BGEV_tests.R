########################
# TESTS
#######################
source("package/bgev_functions.R")
source("package/distCheck.R")
library(fBasics)



# MAKE A TEST SUPPORT FUNCTION AND LOOK AT THE VALUES OF THE DENSITY, MAKE A PLOT, ETC.
# NEGATIVE VALUES IN THE VARIANCE IS BAD, WHEN USING THE INTEGRATION.

#-------------------------------------------
# dbgev
# mu in R; sigma > 0; xi in R ;  delta > -1;
for(i in 1:100){
mu = runif(1,-2,2)
sigma = runif(1,0.1,2)
#xi = 0.1
xi = runif(1,0.3,0.9) * sign(runif(1,-1,1))
delta = 1#runif(1,-0.6,2)
support = bgev.support(mu, sigma, xi, delta)
var.exists = ( xi != 0) & (xi < (delta + 1)/2)
#support = c(-50,50)
ret = distCheck(fun="bgev", n = 2000, 
          support.lower = support[1], support.upper = support[2], subdivisions = 5000,
          mu = mu, sigma = sigma, xi = xi, delta = delta, var.exists = var.exists, print.result = TRUE)
if(ret$test1.density$error.check == FALSE)
  break
print(rbind(unlist(ret$test3.mean.var$computed),c(ret$test3.mean.var$expected$mean$value,ret$test3.mean.var$expected$var,ret$test3.mean.var$expected$log)))
print(bgev.mean(mu = mu, sigma = sigma, xi = xi, delta = delta))
}


mean(rbgev(n = 100000, mu = mu, sigma = sigma, xi = xi, delta = delta))




























# LIMIT TESTS 
# y in R; mu in R; sigma > 0; xi in R ;  delta > -1;
size = 1e5
limitvalues = matrix(c(-size, size,
                       -size, size,
                       size, size,
                       -size, size,
                        0, size), 5, 2, byrow = TRUE)
rownames(limitvalues) = c("y","mu","sigma","xi","delta")
colnames(limitvalues) = c("liminf","limsup")
for(y in limitvalues[1,]){
  for(mu in limitvalues[2,]){
    for(sigma in limitvalues[3,]){
      for(xi in limitvalues[4,]){
        for(delta in limitvalues[5,]){
          print(paste("(y,mu,sigma,xi,delta) = (",y,mu,sigma,xi,delta,") and dbgev = ",try(dbgev(y, mu, sigma, xi, delta)), sep = " " ))
        }
      }
    }
  }
}


dbgev(0, mu = 0, sigma = 0.0001, xi = 0, delta = -0.9) # INVESTIGAR PQ TEM ESSE NA AQUI... ACHO 
# QUE TEM UM INTERTRAVAMENTO ENTRE Y e Delta !!

dbgev(0, mu = 0, sigma = 0.0001, xi = 0, delta = 0)

dbgev(-1, mu = 0, sigma = 0.0001, xi = 0, delta = -0.999)

dbgev(Inf, mu = 0, sigma = 0.0001, xi = 1e5, delta = 0)




# qbgevd
# y in R; mu in R; sigma > 0; xi in R ;  delta > -1;
size = 1e5
limitvalues = matrix(c(1/size, 1-1/size,
                       -size, size,
                       size, size,
                       -size, size,
                       -1 + 1/size, size), 5, 2, byrow = TRUE)
rownames(limitvalues) = c("p","mu","sigma","xi","delta")
colnames(limitvalues) = c("liminf","limsup")
for(p in limitvalues[1,]){
  for(mu in limitvalues[2,]){
    for(sigma in limitvalues[3,]){
      for(xi in limitvalues[4,]){
        for(delta in limitvalues[5,]){
          print(paste("(p,mu,sigma,xi,delta) = (",p,mu,sigma,xi,delta,") and qbgevd = ",try(qbgevd(p, mu, sigma, xi, delta)), sep = " " ))
        }
      }
    }
  }
}

