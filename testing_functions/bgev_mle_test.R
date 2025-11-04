source("cran_versions/bgev_version_0.2/bgev/R/bgev_functions.R")
source("cran_versions/bgev_version_0.2/bgev/R/dist_check.R")
library(fBasics)
library(DEoptim)

####
# TESTE COM PARAMETROS FIGURA QUALIFICACAO DO MATHEUS. ESTUDO SIMULACAO. 
mc = 4
n = 100
parsets = expand.grid(c(-1,0,1),c(1,2),c(-1,1),c(0,1))
colnames(parsets) = c("mu","sigma","xi","delta")
resmc = array(NA,dim = c(mc,4,nrow(parsets)))
for(i in 1:nrow(parsets)){
  print(i)
  par = parsets[i,]
  for(j in 1:mc){
    x = rbgev(n, mu = par$mu,sigma = par$sigma,xi = par$xi,delta = par$delta)
    res = try(bgev.mle(x, method_envstats = "pwme"), silent = TRUE)
    
    # repeat if fail
    while(inherits(res, "try-error")){
      x = rbgev(n, mu = par$mu,sigma = par$sigma,xi = par$xi,delta = par$delta)
      res = try(bgev.mle(x, method_envstats = "pwme"), silent = TRUE)
    }
        resmc[j,,i] = res$par
  } 
}


resmc.mean = t(apply(resmc, MARGIN=c(2,3), FUN=mean))
resmc.sd = t(apply(resmc, MARGIN=c(2,3), FUN=sd))


for(i in 1:nrow(parsets)){
  mc.report = cbind(unlist(parsets[i,]),resmc.mean[i,],resmc.sd[i,])
  colnames(mc.report) = c("True", "Mean", "Sd")
  print(round(mc.report,3))
}
