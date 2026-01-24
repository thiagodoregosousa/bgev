test_that("If bgev(mu, sigma, xi, delta = 0) 
EnvStats::dgevd(location = mu, scale = sigma, shape = -xi)", {

   x_values = seq(-5,5,0.1)
   mu = runif(1,-2,2)
   sigma = runif(1,0.1,2)
   xi = runif(1,0.1,3)
   density_bgev = dbgev(x_values, mu = mu, sigma = sigma, xi = xi, delta = 0)
   density_gev = EnvStats::dgevd(x_values, location = mu, scale = sigma, shape = -xi)
   for(i in 1:length(density_bgev)){
     expect_equal(density_bgev[i], density_gev[i])
   }
   expect_equal(density_bgev, density_gev)
})