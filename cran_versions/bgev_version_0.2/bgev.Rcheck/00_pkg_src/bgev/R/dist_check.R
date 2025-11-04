distCheck = function (fun = "norm", n = 1000, robust = FALSE, subdivisions = 1500, 
                      support.lower = -Inf, support.upper = Inf, var.exists = TRUE, print.result = TRUE,
                      ...) 
{
  # ============= 
  # FUNCTION:
  # ============= 
  # test implementation of continuous distributions in R using its naming convention
  # Disclaimer: this function is just an adaptation of the R function fBasics::distCheck
  # ============= 
  # ARGUMENTS:
  # ============= 
  # fun: name of the function to be tested, e.g., "norm", "gev", "exp"
  # n: size of the sample when generating random values through "rfun"
  # robust: for computing mean and variance in a robust way when evaluating mean and variance
  # subdivisions: used for numerical integration when using "dfun"
  # support.lower and support.upper: lower and upper limit of the support of the distribution
  # var.exists: Boolean indicating if variance exists (useful for gev, bimodal gev, stable etc)
  # ============= 
  # RETURN:
  # ============= 
  # A list containing values computed, expected and the report of checks comparing
  # computed and expected values 
  # ============= 
  # DESCRIPTION: 
  # ============= 
  # to allow for custom distribution support and returning the test values 
  # into an object for future use
  # Description:
  # Tests density, quantile, probability and random function (continuous distributions)
  # 3 tests are done: 
  # -------------
  # test1.density
  # -------------  
  # test if the density function integrates out to 1. Notice here
  # that it would be advisable giving better limiting lower and upper bounds
  # for integration in case the function is not defined in some region, e.g., GEV or Bimodal GEV
  # Bimodal GEV distribution.
  # -------------
  # test2.quantile.cdf
  # -------------
  # Sample moments using the RNG of the variable and comparing with 
  # the value obtained from integrating the density. 
  # CAREFULL FOR DISTRIBUTIONS WHICH DOES NOT HAVE WELL DEFINED FIRST AND/OR SECOND MOMENTS.
  #   -------------
  #   test3.mean.var
  #   -------------
  # compute mean and variance using integral of density and using "rfun" 
  # and compare both at the end.
  # ============= 
  
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
  ret$test1.density$computed = integrate(dfun, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)
  ret$test1.density$expected = 1
  ret$test1.density$error.check = (abs(ret$test1.density$computed[[1]] - 1) < 0.01)
  
  # test2.quantile.cdf
  p = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  ret$test2.quantile.cdf$computed = pfun(qfun(p, ...), ...)
  ret$test2.quantile.cdf$expected = p  
  RMSE = sd(ret$test2.quantile.cdf$computed - ret$test2.quantile.cdf$expected)
  ret$test2.quantile.cdf$error.check = (abs(RMSE) < 1e-04)
  
  # test3.mean.var  
  # computed using "rfun"
  r = rfun(n = n, ...)
  if (!robust) {
    sample.mean = mean(r)
    sample.var = var(r)
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
  exact.mean = integrate(fun1, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)
  exact.second.moment = integrate(fun2, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                  stop.on.error = FALSE, ...)
  exact.log.moment = integrate(fun3, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
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