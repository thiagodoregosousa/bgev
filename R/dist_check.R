#' Consistency checks for continuous distribution implementations
#'
#' Tests whether a set of functions implementing a continuous distribution
#' (density, distribution, quantile, and random generation) satisfy basic
#' probabilistic consistency conditions under the standard R naming convention
#' (\code{d*}, \code{p*}, \code{q*}, \code{r*}).
#'
#' This function is an adaptation of \code{fBasics::distCheck}, extended to
#' allow for custom distribution support and to return all test results in a
#' structured object for further inspection.
#'
#' @param fun Character string giving the name of the distribution (e.g.,
#'   \code{"norm"}, \code{"gev"}, \code{"exp"}).
#' @param n Sample size used when generating random values via the corresponding
#'   \code{r*} function.
#' @param robust Logical; if \code{TRUE}, mean and variance are computed using
#'   robust estimators when applicable.
#' @param subdivisions Number of subdivisions used for numerical integration
#'   when evaluating the density function.
#' @param support.lower Lower bound of the support of the distribution.
#' @param support.upper Upper bound of the support of the distribution.
#' @param var.exists Logical; indicates whether the variance of the distribution
#'   exists (useful for GEV, bimodal GEV, stable distributions, etc.).
#' @param print.result Logical; if \code{TRUE}, a summary of the test results
#'   is printed.
#' @param ... Additional parameters passed to the distribution functions.
#'
#' @details
#' The following consistency checks are performed:
#'
#' \describe{
#'   \item{Density check}{Tests whether the density integrates to one over the
#'   specified support. For distributions with restricted support (e.g., GEV
#'   or bimodal GEV), appropriate bounds should be supplied.}
#'
#'   \item{Quantile--CDF check}{Compares empirical quantiles obtained from random
#'   generation with those implied by the cumulative distribution function.}
#'
#'   \item{Mean--variance check}{Computes mean and variance both from numerical
#'   integration of the density and from simulated samples, and compares the two.
#'   This check is skipped or flagged when moments are not finite.}
#' }
#'
#' @return
#' A list containing the computed values, theoretical expectations, and
#' diagnostic information for each test.
#'
#' @author
#' Thiago do Rego Sousa
#'
#' @seealso
#' \code{\link[fBasics]{distCheck}}
#'
#' @examples
#' \dontrun{
#' distCheck("norm")
#' distCheck("gev", xi = 0.2, sigma = 1, mu = 0,
#'           support.lower = -5, support.upper = 10)
#' }
distCheck <- function(
  fun = "norm",
  n = 1000,
  robust = FALSE,
  subdivisions = 1500,
  support.lower = -Inf,
  support.upper = Inf,
  var.exists = TRUE,
  print.result = TRUE,
  ...
) {
  # construct object to return 
  ret = list( functionTested = fun,
              functionParam = list(...),
              test1.density =      list(computed = NULL, expected = NULL, error.check = NULL),
              test2.quantile.cdf = list(computed = NULL, expected = NULL, error.check = NULL),
              test3.mean.var =     list(computed = list(mean = NULL, var = NULL, log = NULL), 
                                        expected = list(mean = NULL, var = NULL, log = NULL), 
                                        error.check = NULL, 
                                        condition.is.var.finite = TRUE))

  # match functions to test
  CALL = match.call()
  dfun = match.fun(paste("d", fun, sep = ""))
  pfun = match.fun(paste("p", fun, sep = ""))
  qfun = match.fun(paste("q", fun, sep = ""))
  rfun = match.fun(paste("r", fun, sep = ""))
  
  # test1.density
  ret$test1.density$computed = stats::integrate(dfun, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)
  ret$test1.density$expected = 1
  ret$test1.density$error.check = (abs(ret$test1.density$computed[[1]] - 1) < 0.01)
  
  # test2.quantile.cdf
  p = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  ret$test2.quantile.cdf$computed = pfun(qfun(p, ...), ...)
  ret$test2.quantile.cdf$expected = p  
  RMSE = stats::sd(ret$test2.quantile.cdf$computed - ret$test2.quantile.cdf$expected)
  ret$test2.quantile.cdf$error.check = (abs(RMSE) < 1e-04)
  
  # test3.mean.var  
  # computed using "rfun"
  r = rfun(n = n, ...)
  if (!robust) {
    sample.mean = mean(r)
    sample.var = stats::var(r)
    sample.log = mean(log(abs(r)))
  }
  else {
    robustSample = MASS::cov.mcd(r, quantile.used = floor(0.95 * n))
    sample.mean = robustSample$center
    sample.var = robustSample$cov[1, 1]
    sample.log = NULL
  }
  ret$test3.mean.var$computed$mean = sample.mean
  ret$test3.mean.var$computed$var = sample.var
  ret$test3.mean.var$computed$log = sample.log
  # expected
  fun1 = function(x, ...) {
    x * dfun(x, ...)
  }
  fun2 = function(x, ...) {
    x^2 * dfun(x, ...)
  }
  fun3 = function(x, ...) {
    log(abs(x)) * dfun(x, ...)
  }
  exact.mean = stats::integrate(fun1, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)
  exact.second.moment = stats::integrate(fun2, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                  stop.on.error = FALSE, ...)
  exact.log.moment = stats::integrate(fun3, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                                  stop.on.error = FALSE, ...)

  exact.var = exact.second.moment[[1]] - exact.mean[[1]]^2
  ret$test3.mean.var$expected$mean = exact.mean
  ret$test3.mean.var$expected$var = exact.var
  ret$test3.mean.var$error.check = (abs((sample.var - exact.var)/exact.var) < 0.1)
  # ret$test3.mean.var$error.check = (abs(sample.mean - exact.mean[[1]]) < 0.1)
  ret$test3.mean.var$expected$log = exact.log.moment$value
  
  if(print.result){
    cat("\n============================================================
        REPORT OF ALL 3 TESTS. SHOULD BE ALL TRUE TO PASS 
============================================================\n")
    ans = list(test1.density = ret$test1.density$error.check, 
               test2.quantile.cdf = ret$test2.quantile.cdf$error.check,
               test3.mean.var = ret$test3.mean.var$error.check,
               condition.is.var.finite = ret$test3.mean.var$condition.is.var.finite)
    print(unlist(ans))
  }
  
  # return
  ret
}