\name{bgev.mle}
\alias{bgev.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation of bimodal GEV distribution based on real data 
that appears to be bimodal.
}
\description{
Finds the maximum likelihood estimators of the bimodal GEV distribution 
parameters by minimizing the log-likelihood function of the data.

The estimation procedure is a divided into three steps:

1) Use egevd for starting values: We estimate the data using a unimodal GEV distribution (function egevd from package EnvStats). This is justified 
by the fact that the bgev distribution with delta = 0 is a unimodal GEV distribution with parametesr mu sigma and xi.  

2) DEOptim for more adequate starting values: Use the starting values to create lower and upper bounds and feed it into DEOptim to get starting values for all
parameters simutaneously.

3) Use optim for local search with starting values for DEOptim.

}
\usage{
bgev.mle(x, method_envstats = "mle", deoptim.itermax = 200, 
optim.method = "L-BFGS-B",start = NULL, lower = NULL, upper = NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a unidimensional vector containing observations to estimate a bimodal GEV 
distribution
}
\item{start}{
a vector of dimension 4 containing the starting values for each of the bimodal GEV
parameters where the search is going to take place. 
}

  \item{lower}{
a vector of dimension 4 containing the lower limit for each of the bimodal GEV
parameters where the search is going to take place. 
}
  \item{upper}{
a vector of dimension 4 containing the upper limit for each of the bimodal GEV
parameters where the search is going to take place. 
}
 \item{method_envstats}{
method to use in the egevd function from envstats package.
}
 \item{deoptim.itermax}{
number of iterations in the deoptim step.
}
  \item{optim.method}{
maximum number of interations when finding a good starting value for the optimization algorithm.
}
}
\value{
A list with components returned by the optim R function for which the 
main values are
\item{par }{best parameters of bimodal gev fitting the data}
\item{value }{value of the log-likelihood corresponding to 'par'}
}
\references{
OTINIANO, Cira EG, et al. (2023). A bimodal model for extremes data. Environmental and Ecological Statistics, 1-28. http://dx.doi.org/10.1007/s10651-023-00566-7
}
\author{
\packageAuthor{bgev}
}
\examples{
# generate 100 values distributed according to a bimodal GEV
x = rbgev(50, mu = 0.2, sigma = 1, xi = 0.5, delta = 0.2) 
# estimate the bimodal GEV parameters using the generated data
bgev.mle(x)
}
\keyword{ bgev }
\keyword{ gev }
\keyword{ bimodal }
\keyword{ mle }
\keyword{ estimation }