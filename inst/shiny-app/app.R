#################################################
# Install once if needed
#################################################
# install.packages("shiny")

library(shiny)
library(bgev)

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
