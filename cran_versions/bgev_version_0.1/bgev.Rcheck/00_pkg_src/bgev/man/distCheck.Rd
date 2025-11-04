\name{distCheck}
\alias{distCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check implementation of a distribution in \R. 
}
\description{
This is just a refactorization of 
the original \R function \code{distCheck} to allow customization for heavy tailed distributions
and distribution with constrained support. Assumes one has implemented density, 
probability, quantile and random generation for a new distribution. 
Originally implemented to be used with the bimodal GEV distribution
}
\usage{
distCheck(fun = "norm", n = 1000, robust = FALSE, subdivisions = 1500, 
support.lower = -Inf, support.upper = Inf, 
var.exists = TRUE, print.result = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{
Name of distribution. E.g, "bgev" for bimodal GEV, "gev" for GEV, "norm" for normal. 
}
  \item{n}{
size of the sample when generating random values through "rfun"
}
  \item{robust}{
 Boolean. Used for computing mean and variance in a robust way when evaluating mean 
 and variance.
}
  \item{subdivisions}{
 used for numerical integration when using "dfun"
}
  \item{support.lower}{
lower limit of the support of the distribution
}
  \item{support.upper}{
upper limit of the support of the distribution
}
  \item{var.exists}{
Boolean indicating if variance exists (useful for gev, bimodal gev, stable etc)
}
  \item{print.result}{
Boolean indicating weather to print results of the tests on screen.
}
  \item{\dots}{
parameters passed to the distribution, e.g., mu and sd for normal, mu, sigma, xi
and delta for bimodal GEV.
}
}
\details{
All tests results are organized into a list to make it easier to find where 
results went wrong and to make it scalable for testing distributions for 
several parameters. 
}
\value{
A list of lists with test results


\code{list( functionTested = fun,
              functionParam = list(...),}\cr
              \code{test1.density =      list(computed = NULL, expected = NULL, error.check = NULL),}\cr
              \code{test2.quantile.cdf = list(computed = NULL, expected = NULL, error.check = NULL),}\cr
              \code{test3.mean.var =     list(computed = list(mean = NULL, var = NULL, log = NULL), }\cr
                                       \code{ expected = list(mean = NULL, var = NULL, log = NULL),}\cr
                                       \code{  error.check = NULL, }\cr
                                       \code{  condition.is.var.finite = TRUE))}
\item{functionParam }{Additional parameters passed to the distribution with the ... argument}
\item{test1.density }{results for the density function test. A list with the 'computed' value 
using the density function and numerical integration, the 'expected' value (which is 1) and 
the result of the 'error.check' comparing both values.}
\item{test2.quantile.cdf }{results for the quantile function test. A list with the 'computed' value 
using the quantile function at chosen probabilites on [0,1], the 'expected' value (which is 1) and 
the result of the 'error.check' comparing both values.}

\item{test3.mean.var }{results for the random generation and the density function. 
A list with the 'computed' values using the random quantities (computes mean(X), variance(X) and E(log(abs(X))) where X follows
the distribution being tested. The 'expected' is a list containing the same quantities computed using numerical integration}

\item{error.check }{check if expected variance and computed variance are close enough.}

\item{condition.is.var.finite }{For some distributions there are known conditions for the variance to be finite. This quantity
is expected to be evaluated outside the DistCheck function and passed through the argument 'var.exists'.}

}
\references{
Otiniano, Cira EG, et al. (2023). A bimodal model for extremes data. Environmental and Ecological Statistics, 1-28. http://dx.doi.org/10.1007/s10651-023-00566-7
}
\author{
\packageAuthor{bgev}
}
\examples{
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
}
\keyword{ distribution }
\keyword{ test }
\keyword{ implementation }