\name{likbgev}
\alias{likbgev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log likelihood function for the bimodal GEV distribution.
}
\description{
Uses the density function to evaluate the likelihood. This is useful for the
'bgev.mle' function.
}
\usage{
likbgev(y, theta = c(1, 1, 0.3, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a unidimensional vector containing the points to compute the log likelihood
}
  \item{theta}{
bimodal GEV parameters as an ordered pair \code{theta = (mu,sigma,xi,delta)}
}
}
\value{
a unidimensional vector containing the computed log likelihood for \code{y}.
}
\references{
Otiniano, Cira EG, et al. (2023). A bimodal model for extremes data. Environmental and Ecological Statistics, 1-28. http://dx.doi.org/10.1007/s10651-023-00566-7
}
\author{
\packageAuthor{bgev}
}
\examples{
# get random points from bimodal GEV
y = rbgev(100, mu = 1, sigma = 1, xi = 0.3, delta = 2)

# compute log-likelihood
likbgev (y, theta = c(1, 1, 0.3, 2)) 
}
\keyword{ bgev }
\keyword{ gev }
\keyword{ bimodal }
\keyword{ likelihood }