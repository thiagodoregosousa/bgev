# Histogram-based heuristic for locating mu in a bimodal sample, extracted from
# R/bgev_diagnostics.R (moved out of R/ because it isn't package API and its
# trailing demo referenced undefined globals at top level, breaking package
# load). Not part of the package; run interactively from the repo root.

source("R/bgev_domain.R")
source("R/bgev_distribution.R")

estimate_mu_from_histogram <- function(x,
                                       binwidth = NULL,
                                       breaks = NULL,
                                       interior_frac = 0.2,   # ignore 20% on each side
                                       min_drop_frac = 0.05,  # ignore tiny drops
                                       plot = FALSE) {

  # ---------------------------
  # 1. Build histogram
  # ---------------------------

  if (!is.null(binwidth)) {
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    breaks <- seq(xmin, xmax + binwidth, by = binwidth)
  }

  if (is.null(breaks)) {
    h <- hist(x, plot = plot)
  } else {
    h <- hist(x, breaks = breaks, plot = plot)
  }

  counts <- h$counts
  mids   <- h$mids
  n_bins <- length(counts)

  if (n_bins < 5) {
    stop("Not enough bins to compute interior drops.")
  }

  # ---------------------------
  # 2. Compute drops
  # ---------------------------

  diffs <- diff(counts)

  # ---------------------------
  # 3. Define interior region
  # ---------------------------

  lower_idx <- floor(n_bins * interior_frac)
  upper_idx <- ceiling(n_bins * (1 - interior_frac))

  # Keep only interior drops
  candidate_idx <- which(diffs < 0 &
                           seq_along(diffs) > lower_idx &
                           seq_along(diffs) < upper_idx)

  if (length(candidate_idx) == 0) {
    return(list(
      mu_estimate = NA,
      message = "No interior drop detected."
    ))
  }

  # ---------------------------
  # 4. Remove very small drops
  # ---------------------------

  max_count <- max(counts)
  min_required_drop <- min_drop_frac * max_count

  candidate_idx <- candidate_idx[abs(diffs[candidate_idx]) >= min_required_drop]

  if (length(candidate_idx) == 0) {
    return(list(
      mu_estimate = NA,
      message = "Interior drops too small (noise-level)."
    ))
  }

  # Choose strongest interior drop
  best_idx <- candidate_idx[which.min(diffs[candidate_idx])]

  mu_estimate <- mids[best_idx]

  return(list(
    mu_estimate = mu_estimate,
    drop_value = diffs[best_idx],
    drop_index = best_idx,
    counts = counts,
    mids = mids,
    breaks = h$breaks,
    binwidth = diff(h$breaks)[1]
  ))
}

# Demo
mu <- 0; sigma <- 1; xi <- 0.5; delta <- 1
x <- rbgev(n = 1000, mu = mu, sigma = sigma, xi = xi, delta = delta)
estimate_mu_from_histogram(x)
