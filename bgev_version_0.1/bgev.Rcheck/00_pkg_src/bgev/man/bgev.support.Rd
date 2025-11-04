\name{bgev.support}
\alias{bgev.support}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Support of the bimodal GEV distribution
}
\description{
When the shape parameter \code{xi} is different from zero, the support
is truncated either at the left or at the right side of the real. 
Considering the support is particularly useful to estimating momoments and 
to compute the likelihood function.
}
\usage{
bgev.support(mu, sigma, xi, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
location parameter
}
  \item{sigma}{
shape parameter
}
  \item{xi}{
shape parameter
}
  \item{delta}{
bimodality parameter
}
}
\value{
It returns a vector representing the interval for which the density function 
is not null. 
}
\references{
Otiniano, Cira EG, et al. (2023). A bimodal model for extremes data. Environmental and Ecological Statistics, 1-28. http://dx.doi.org/10.1007/s10651-023-00566-7
}
\author{
\packageAuthor{bgev}
}
\note{
The Support of the bimodal GEV distribution is given by the following equation:

\deqn{
    \begin{cases}
     x \in \mathbb {R} : x \geq \mathbf{sng} \left(-\frac{\sigma}{\xi}\right) | \frac{\sigma}{\xi}| ^{\frac{1}{\delta+1}}+\mu, &  \xi \neq 0 \\
     x \in \mathbb {R}, & \xi =0.
    \end{cases}
}


}
\examples{
# Computes the support of a specific bimodal GEV distribution
support = bgev.support(mu=1, sigma=10, xi=0.3, delta=2)
}
\keyword{ bgev }
\keyword{ gev }
\keyword{ bimodal }
\keyword{ support }