\name{bgev-package}
\alias{bgev-package}
\alias{bgev}
\alias{dbgev}
\alias{pbgev}
\alias{qbgev}
\alias{rbgev}
\docType{package}
\title{
\packageTitle{bgev}
}
\description{
\packageDescription{bgev}
}
\details{
Density, distribution function, quantile function and random generation of 
bimodal GEV distribution with location parameter. In addition, maximum likelihood
estimation based on real data is also provided.
}
\usage{
dbgev(y, mu = 1, sigma = 1, xi = 0.3, delta = 2)
pbgev(y, mu = 1, sigma = 1, xi = 0.3, delta = 2)
qbgev(p, mu = 1, sigma = 1, xi = 0.3, delta = 2)
rbgev(n, mu = 1, sigma = 1, xi = 0.3, delta = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a unidimensional vector containing the points to compute the density (dbgev) or  
the probability (pbgev) froma bimodal GEV distribution with parameters
\code{mu}, \code{sigma}, \code{xi} and \code{delta}.
}
  \item{p}{
a unidimensional vector containing the probabilities 
used to compute the quantiles
}
  \item{n}{
an integer describing the number of observations to generate random bimodal
GEV observations
}
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
\code{dbgev} gives the density, \code{pbgev} gives the distribution function, 
\code{qbgev} gives the quantile function, and \code{rbgev} generates random
bimodal GEV observations. 
}
\author{
\packageAuthor{bgev}

Maintainer: \packageMaintainer{bgev}
}
\note{

The probability density of a GEV random variable; \eqn{Y \sim F_{\xi, \sigma, \mu}} is given by: 

\deqn{
	f_{\xi, \mu, \sigma}(y)=
		\begin{cases}
		\dfrac{1}{\sigma} \left[ 1+ \xi \left(\dfrac{y-\mu}{\sigma}\right) \right]^{(-1/\xi) -1} \exp\left\{- \left[1+\xi\left(\dfrac{y-\mu}{\sigma}\right)\right]^{-1/\xi}\right\} ,& \text{if } \xi \ne 0 \\
		\dfrac{1}{\sigma} \exp \left\{ - \left(  \dfrac{y-\mu}{\sigma}\right)  - \exp \left[ - \left(  \dfrac{y-\mu}{\sigma}\right)  \right]  \right\},              & \text{if }  \xi = 0 ,
	\end{cases}
}
where \eqn{\xi} and \eqn{\sigma} are the shape parameters and \eqn{\mu} is the location parameter.

The bimodal Generalized Extreme Value (GEV) model, denoted as BGEV, consists of composing the distribution of a random variable following the GEV distribution with a location parameter \eqn{\mu=0}, i.e., \eqn{Y \sim F_{\xi, 0, \sigma}}, with the transformation \eqn{T_{\mu, \delta}} defined below. Thus, the cumulative distribution function (CDF) of a random variable BGEV, denoted as \eqn{X \sim F_{BG_{\xi, \mu, \sigma, \delta}}}, is given by:

\deqn{F_{BG_{\xi,\mu,\sigma, \delta}}(x) =  F_{\xi, 0, \sigma}(T_{\mu, \delta}(x)),}

where the function \eqn{T_{\mu, \delta}} is defined as:

\deqn{T_{\mu, \delta}(x)=\left(  x - \mu  \right) \left| x -\mu \right| ^{\delta}, \quad \delta > -1, \quad \mu \in \mathbb{R}.}

Moreover, the function \eqn{T} is invertible, with the inverse given by:

\deqn{T^{-1}_{\mu, \delta}(x) = \text{sng}(x) |x|^{\dfrac{1}{\left( \delta +1 \right) }} + \mu.}

Additionally, it is differentiable, and its derivative has the following form:

\deqn{T'_{\mu, \delta}(x) = (\delta + 1 ) |x - \mu|^{\delta}.}


Its probability density function \eqn{X\sim F_{BG_{\xi,\mu,\sigma, \delta}}} is given by

\deqn{
	f_{BG_{\xi,\mu,\sigma, \delta}} (x)= \begin{cases}
		\dfrac{1}{\sigma} \left[ 1+ \xi \left(\dfrac{T_{\mu, \delta}(x)}{\sigma}\right) \right]^{(-1/\xi) -1} \exp\left[- \left[1+\xi\left(\dfrac{T_{\mu, \delta}(x)}        {\sigma}\right)\right]^{-1/\xi}\right] T'_{\mu, \delta}(x)
		, & \xi \neq 0 \\
		\dfrac{1}{\sigma} \exp \left( - \dfrac{T_{\mu, \delta}(x)}{\sigma} \right) \exp \left[- \exp \left( - \dfrac{T_{\mu, \delta}(x)}{\sigma}\right)  \right] T'_{\mu, \delta}(x), &  \xi=0.
	\end{cases}
}

}
\references{
Otiniano, Cira EG, et al. (2023). A bimodal model for extremes data. Environmental and Ecological Statistics, 1-28. http://dx.doi.org/10.1007/s10651-023-00566-7
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
\keyword{ distribution }
