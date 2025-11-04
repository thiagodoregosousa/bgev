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
}
\usage{
bgev.mle(x, lower = c(-3, 0.1, -3, -0.9), upper = c(3, 3, 3,3), itermax = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a unidimensional vector containing observations to estimate a bimodal GEV 
distribution
}
  \item{lower}{
a vector of dimension 4 containing the lower limit for each of the bimodal GEV
parameters where the search is going to take place. 
}
  \item{upper}{
a vector of dimension 4 containing the upper limit for each of the bimodal GEV
parameters where the search is going to take place. 
}
  \item{itermax}{
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