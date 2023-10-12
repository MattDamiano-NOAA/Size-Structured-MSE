# Dolphinfish MSE
# OM wrapper for Shiny app development
# Last date modified: 10-11-2023 by M. Damiano

# set your OM directory to where you store all three OM files

om.dir = "C:/Users/matt.damiano/Desktop/Dolphinfish MSE/OM files"

# set working directory

setwd(om.dir)

# run the Functions file - only need to do this once
source("Functions.r")

# set iteration dimension
n.iter = 100

# create an array to store your spawning stock biomass "estimates"
ssb.array <- array(NA, dim = c(37, n.iter)) # since we haven't called the input file yet, 
# and which contains the other model dimensions, let's give it the length of 
# the annual time series (37), and the number of iterations we want to run the 
#simulation for (100/'n.iter')

# stochastic recruitment? (TRUE/FALSE)
rec.sto = TRUE # try this with and without stochasticity
# without it, all simulations should be deterministic (same result)

iter = 1

while(iter<=n.iter){
  
  source("input.data.r") # run input file
  source("data.gen.R") # produces results for OM
  
  ssb.array[,iter] = PopDy$SSB # store SSB for each iteration
  
  iter = iter+1 # update counter 
}

# jitter plot to visualize (stochastic) results
windows(height=10, width = 10)
plot(seq(1986,2022,1),rowMeans(ssb.array, dim=1), type="l",ylim=c(0,2e9), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="SSB", xlab="Year")
for (i in 1:n.iter){
  points(seq(1986,2022,1),jitter(ssb.array[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(ssb.array[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
}