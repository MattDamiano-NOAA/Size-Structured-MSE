# Dolphinfish MSE
# OM wrapper for Shiny app development
# Last date modified: 10-11-2023 by M. Damiano

om.dir = "./GitHub/Size-Structured-MSE"
setwd(om.dir)
source("Functions.r") # run the Functions file - only need to do this once

# Main display calculations
n.iterations = 5
ssb.array <- array(NA, dim = c(37, n.iterations))
abundance.array <- array(NA, dim = c(37, 20, 4, 7, n.iterations))
rec.sto = TRUE # stochastic recruitment? (TRUE/FALSE)
iteration = 1

while(iteration<=n.iterations){
  source("input.data.r") # run input file
  source("data.gen.R") # produces results for OM
  ssb.array[,iteration] = PopDy$SSB # store SSB for each iteration
  abundance.array[,,,,iteration] = PopDy$Abundance
  iteration = iteration + 1 # update counter
}

dev.new(width=50, height=20, unit="in") # open new window

plot(
  seq(1986, 2022, 1),
  rowMeans(ssb.array, dim = 1),
  type = "l",
  ylim = c(min(ssb.array), max(ssb.array)),
  col = rgb(160, 32, 240, 255, maxColorValue=255),
  lwd = 3,
  yla = "SSB",
  xlab = "Year"
)
for (i in 1:n.iterations){
  points(
    seq(1986, 2022, 1),
    jitter(ssb.array[,i], factor = 2),
    type = "l",
    col = rgb(
      160, 32, 240,
      (255 - abs(round(dim(ssb.array[,1:n.iterations])[2] / 2) - i))/9,
      maxColorValue = 255
    ),
    lwd = 2)
}
