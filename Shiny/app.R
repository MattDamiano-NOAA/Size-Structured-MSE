#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(tidyr)
source("R/functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Dolphinfish MSE with Modules"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("iterations",
                  "Number of iterations:",
                  min = config::get('min_iterations'),
                  max = config::get('max_iterations'),
                  value = config::get('iterations')),
      checkboxInput("stochasticity",
                    "Stochasticity?",
                    value = TRUE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("ssbPlot"),
      plotOutput("abundancePlot"),
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  initial_year <- 1986
  terminal_year <- 2022
  output$ssbPlot <- renderPlot({
    # Input parameters
    n.iterations <- input$iterations

    # Main display calculations
    ssb.array <- array(NA, dim = c(37, n.iterations))
    iteration_cols <- round(dim(ssb.array[,1:n.iterations])[2] / 2)
    iteration = 1

    rec.sto = input$stochasticity # stochastic mean recruitment?
    while (iteration <= n.iterations) {
      source("R/data.params.year_size.R")
      source("R/data.params.selectivity.R")
      source("R/data.params.movement.R")
      source("R/data.params.weight_length_relationship.R")
      source("R/data.params.maturity.R")
      source("R/data.params.natural_mortality.R")
      source("R/data.params.fishing_mortality.R")
      source("R/data.params.growth.R")
      source("R/data.params.recruitment.R")
      source("R/data.params.pop_dynamics.R")
      source("R/data.params.cpue.R")
      source("R/data.params.survey.R")
      source("R/data.gen.population_dynamics.R") # produces results for OM
      ssb.array[,iteration] = PopDy$SSB # store SSB for each iteration
      iteration = iteration + 1 # update counter
    }

    plot(
      seq(initial_year, terminal_year, 1),
      rowMeans(ssb.array, dims = 1),
      type = "l",
      ylim = c(min(ssb.array), max(ssb.array)),
      col = rgb(160, 32, 240, 255, maxColorValue = 255),
      lwd = 3,
      ylab = "SSB",
      xlab = "Year",
      main = "Spawning Stock Biomass"
    )
  })

  output$abundancePlot <- renderPlot({
    # Input parameters
    n.iterations <- input$iterations

    # Main display calculations
    ssb.array <- array(NA, dim = c(37, n.iterations))
    iteration_cols <- round(dim(ssb.array[,1:n.iterations])[2] / 2)
    abundance.array <- array(NA, dim = c(37, 20, 4, 7, n.iterations))
    iteration = 1

    rec.sto = input$stochasticity # stochastic mean recruitment?
    while (iteration <= n.iterations) {
      source("R/data.params.year_size.R")
      source("R/data.params.selectivity.R")
      source("R/data.params.movement.R")
      source("R/data.params.weight_length_relationship.R")
      source("R/data.params.maturity.R")
      source("R/data.params.natural_mortality.R")
      source("R/data.params.fishing_mortality.R")
      source("R/data.params.growth.R")
      source("R/data.params.recruitment.R")
      source("R/data.params.pop_dynamics.R")
      source("R/data.params.cpue.R")
      source("R/data.params.survey.R")
      source("R/data.gen.population_dynamics.R") # produces results for OM
      abundance.array[,,,,iteration] = PopDy$Abundance
      iteration = iteration + 1 # update counter
    }

    plot(
      seq(initial_year, terminal_year, 1),
      rowMeans(abundance.array, dims = 1),
      type = "l",
      ylim = c(min(abundance.array), max(abundance.array)),
      col = rgb(160, 32, 240, 255, maxColorValue = 255),
      lwd = 3,
      ylab = "Abundance",
      xlab = "Year",
      main = "Abundance"
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
