- incluir qqplot para bgev, util para estimacao. 
- put newline at DESCRIPTION with Encoding: UTF-8
- delta upper bound use 20 instead of 2. 
- change start of delta from 0.1 to something using the data
- is it actualy acceptoing lower and upper ?
- gives error
library(bgev)
x = bgev::rbgev(1000, mu = 0, sigma = 1, xi = 0, delta = 9)
bgev::bgev.mle(x, lower = c(-20, 0.001, -20, -0.99), upper = c(20, 10, 20, 20))
- return(1e99) may be bad for some parameters, better to return other values like 1e6
- bgev.mle get new version used in mbev code. 
- bgev.mle currently fails if cannot get starting values for mle. Use my suggestion error message to actually 
do that during estimation, maybe just print a message. 