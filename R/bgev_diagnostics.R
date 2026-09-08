#################################################
# Install once if needed
#################################################
# install.packages("shiny")

library(shiny)

#################################################
# UI
#################################################

ui <- fluidPage(
  
  titlePanel("Interactive BGEV Density Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput("mu", "mu:",
                  min = -20, max = 20, value = 0, step = 0.1),
      
      sliderInput("sigma", "sigma:",
                  min = 0.1, max = 10, value = 1, step = 0.1),
      
      sliderInput("xi", "xi:",
                  min = -5, max = 5, value = 0.5, step = 0.1),
      
      sliderInput("delta", "delta:",
                  min = -2, max = 5, value = 1, step = 0.1),
      
      numericInput("n", "Sample size:",
                   value = 200, min = 10, max = 5000),
      
      actionButton("simulate", "Simulate Sample")
    ),
    
    mainPanel(
      plotOutput("densityPlot")
    )
  )
)

#################################################
# SERVER
#################################################

server <- function(input, output) {
  
  # Reactive dataset
  data_sample <- eventReactive(input$simulate, {
    rbgev(n = input$n,
          mu = input$mu,
          sigma = input$sigma,
          xi = input$xi,
          delta = input$delta)
  }, ignoreNULL = FALSE)
  
  output$densityPlot <- renderPlot({
    
    x <- seq(min(data_sample()) - 2,
             max(data_sample()) + 2,
             length.out = 1000)
    
    y <- dbgev(x,
               mu = input$mu,
               sigma = input$sigma,
               xi = input$xi,
               delta = input$delta)
    
    # Plot histogram
    hist(data_sample(),
         probability = TRUE,
         breaks = 30,
         col = "lightgray",
         border = "white",
         main = "BGEV Density vs Simulated Data",
         xlab = "x")
    
    # Overlay theoretical density
    lines(x, y, col = "red", lwd = 3)
  })
}

#################################################
# RUN APP
#################################################

shinyApp(ui = ui, server = server)













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





# Simulate data with true parameters
x <- rbgev(n = 1000, mu = mu, sigma = sigma, xi = xi, delta = delta)
estimate_mu_from_histogram(x)
