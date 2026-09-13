# Exploratory/dev scratch code for R/bgev_start_params.R, extracted from that file
# because top-level executable statements in R/ run at package load time and were
# breaking devtools::document()/check()/load_all(). Not part of the package; run
# interactively from the repo root.

source("R/bgev_domain.R")
source("R/bgev_distribution.R")
source("R/bgev_estimation.R")
source("R/bgev_start_params.R")

# Quantile identity check: quantile at y = exp(-1) equals mu for bgev
mu = 5 # runif(1,-1000,1000)
sigma =  runif(1,0.01,10)
xi = 0 # runif(1,-10,10)
delta = runif(1,-0.99,10)
qbgev(p = exp(-1), mu = mu, sigma = sigma, xi = xi, delta = delta) == mu

q1 = qbgev(p = exp(-exp(2)), mu = mu, sigma = sigma, xi = xi, delta = delta)
q2 = qbgev(p = exp(-exp(1)), mu = mu, sigma = sigma, xi = xi, delta = delta)
round(1/log(q1/q2, 2) - 1 - delta, 8) == 0
round(sigma == (-q2)^(delta+1), 8) == 0

# Monte Carlo check of mc_estimator's bias/sd at a fixed parameter set
mc_estimator(100000, R = 100, mu = 3, sigma = 2, delta = 2, xi = 0)

# Sanity plot of 2^(-xi)
xi = seq(-3,3,0.01)
plot(xi, 2^(-xi), type = 'l')
abline(h = 1)

# Coverage check: how often bgev_start_region's box contains the true parameters
total <- 100
count <- 0

for (i in 1:total) {

  # Generate true parameters
  mu    <- runif(1, -1000, 100)
  sigma <- runif(1, 0.01, 1000)
  xi    <- runif(1, -10, 10)
  delta <- runif(1, -0.99, 10)

  # Simulate data with true parameters
  x <- rbgev(n = 1000, mu = mu, sigma = sigma, xi = xi, delta = delta)

  # Get estimated starting region
  start_result <- bgev_start_region(x)

  lw <- start_result$lower
  up <- start_result$upper

  # Check if true parameters lie inside bounds
  inside = FALSE
  inside <- (
    lw[1] <= mu    && mu      <= up[1] &&
      lw[2] <= sigma && sigma <= up[2] &&
      lw[3] <= xi    && xi    <= up[3] &&
      lw[4] <= delta && delta <= up[4]
  )

  if(inside == TRUE)
    print(c(i,count))

  if (inside) {
    count <- count + 1
  }
}

# Percentage of times bounds contain true parameters
percentage <- 100 * count / total

cat("Coverage percentage:", percentage, "%\n")
