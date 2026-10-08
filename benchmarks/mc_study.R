# Monte Carlo study for bgev_mle (new API). Purpose: validate the estimator and
# characterise its diagnostics across the regularity regimes, to pave the way to
# CRAN. Simulate from known (mu, sigma, xi, delta), fit, and measure recovery,
# convergence/admissibility rates, and Wald-coverage degradation near the
# parameter-dependent support boundary.
#
# Grid rationale:
#  - xi drives the support non-regularity: xi=0 is the regular control (support
#    = R); |xi|>0 turns on a finite endpoint (sign = which side).
#  - delta>0 only (estimation restricts delta>0): 0.25 ~ near-GEV, 1 moderate
#    bimodal, 3 strong bimodal / heavy tail.
#  - n small->moderate, where small-sample instability is expected.

suppressMessages(devtools::load_all("/Users/thiago/Documents/bgev_github"))
suppressMessages(library(parallel))

mu_true <- 0; sigma_true <- 1
xi_grid    <- c(-0.4, -0.2, 0, 0.2, 0.4)
delta_grid <- c(0.25, 1, 3)
n_grid     <- c(100, 250, 500)
R          <- as.integer(Sys.getenv("MC_R", "300"))  # replications per cell
n_starts   <- 15
z          <- qnorm(0.975) # 95% Wald

out_dir <- "/Users/thiago/Documents/bgev_github/benchmarks/mc_study_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- expand.grid(xi = xi_grid, delta = delta_grid, n = n_grid,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
tasks <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g)
  data.frame(cell = g, grid[g, ], rep = seq_len(R))))

truth_of <- function(xi, delta) c(mu = mu_true, sigma = sigma_true, xi = xi, delta = delta)

fit_one <- function(i) {
  tk <- tasks[i, ]
  truth <- truth_of(tk$xi, tk$delta)
  set.seed(1e6 * tk$cell + tk$rep)
  x <- rbgev(tk$n, mu_true, sigma_true, tk$xi, tk$delta)

  res <- tryCatch(suppressWarnings(bgev_mle(x, n_starts = n_starts)), error = function(e) NULL)
  base <- data.frame(cell = tk$cell, xi = tk$xi, delta = tk$delta, n = tk$n, rep = tk$rep)
  if (is.null(res)) return(cbind(base, failed = TRUE, t(setNames(rep(NA, 14),
    c("mu","sigma","xi_hat","delta_hat","se_mu","se_sigma","se_xi","se_delta",
      "cov_mu","cov_sigma","cov_xi","cov_delta","conv","admissible")))))

  est <- res$par
  # observed information = Hessian of the (continuous) negative log-lik at mle
  se <- rep(NA_real_, 4)
  H <- tryCatch(numDeriv::hessian(function(p) bgev_negative_log_likelihood(x, p), est),
                error = function(e) NULL)
  if (!is.null(H)) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V)) { d <- diag(V); d[d < 0] <- NA; se <- sqrt(d) }
  }
  lo <- est - z * se; hi <- est + z * se
  cov <- as.integer(truth >= lo & truth <= hi)

  cbind(base, failed = FALSE,
        mu = est[1], sigma = est[2], xi_hat = est[3], delta_hat = est[4],
        se_mu = se[1], se_sigma = se[2], se_xi = se[3], se_delta = se[4],
        cov_mu = cov[1], cov_sigma = cov[2], cov_xi = cov[3], cov_delta = cov[4],
        conv = res$convergence, admissible = res$admissible)
}

cat(sprintf("MC: %d cells x %d reps = %d fits on %d cores\n",
            nrow(grid), R, nrow(tasks), max(1, detectCores() - 1)))
t0 <- Sys.time()
rows <- mclapply(seq_len(nrow(tasks)), fit_one, mc.cores = max(1, detectCores() - 1))
raw <- do.call(rbind, rows)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("done in %.1f min\n", elapsed))

saveRDS(raw, file.path(out_dir, "mc_raw.rds"))

# ---- per-cell aggregation ----
agg <- do.call(rbind, by(raw, raw$cell, function(d) {
  ok <- d[d$failed == FALSE, ]
  tr <- truth_of(ok$xi[1], ok$delta[1])
  bias <- c(mean(ok$mu) - tr[1], mean(ok$sigma) - tr[2],
            mean(ok$xi_hat) - tr[3], mean(ok$delta_hat) - tr[4])
  rmse <- c(sqrt(mean((ok$mu - tr[1])^2)), sqrt(mean((ok$sigma - tr[2])^2)),
            sqrt(mean((ok$xi_hat - tr[3])^2)), sqrt(mean((ok$delta_hat - tr[4])^2)))
  data.frame(xi = d$xi[1], delta = d$delta[1], n = d$n[1],
    fail_rate = mean(d$failed), conv_rate = mean(ok$conv == 0, na.rm = TRUE),
    admissible_rate = mean(ok$admissible, na.rm = TRUE),
    bias_mu = bias[1], bias_sigma = bias[2], bias_xi = bias[3], bias_delta = bias[4],
    rmse_mu = rmse[1], rmse_sigma = rmse[2], rmse_xi = rmse[3], rmse_delta = rmse[4],
    cov_mu = mean(ok$cov_mu, na.rm = TRUE), cov_sigma = mean(ok$cov_sigma, na.rm = TRUE),
    cov_xi = mean(ok$cov_xi, na.rm = TRUE), cov_delta = mean(ok$cov_delta, na.rm = TRUE))
}))
rownames(agg) <- NULL

saveRDS(agg, file.path(out_dir, "mc_summary.rds"))
write.csv(agg, file.path(out_dir, "mc_summary.csv"), row.names = FALSE)
cat("\n===== PER-CELL SUMMARY =====\n")
print(round(agg, 3))
cat(sprintf("\nsaved to %s\n", out_dir))
