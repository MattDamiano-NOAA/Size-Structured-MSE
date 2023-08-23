# Size-structured MSE 
# MSE projection file
# M. Damiano
# Last updated: 8-23-2023

# load libraries for running AD Model Builder code
require(PBSadmb)
require(psychTools) 

# Set up your Directories
AM.code.dir = "C:/Users/matt.damiano/Desktop/MARFIN MSE/Assessment model and figure code"
OM.code.dir = "C:/Users/matt.damiano/Desktop/MARFIN MSE/Operating Model and Projection code"
result.dir = "C:/Users/matt.damiano/Desktop/MARFIN MSE/Results"

source('C:/Users/matt.damiano/Desktop/MARFIN MSE/Assessment model and figure code/reptoRlist.r')

# condition the population dynamics - OM set-up
setwd(OM.code.dir)
source("Functions.r")
source("input.data.r")
source("set.pars.r")
source("mse.proj.r")

########################################## simulation paths setting #######################################
scenario         = c('Uber test 6-28-2023') 

sce.dir          = paste(result.dir,'/',scenario,sep='')
input.dir        = paste(result.dir,'/',scenario,'/InputFiles',sep='')
input.season.dir = paste(input.dir,'/Season',sep='')
input.year.dir   = paste(input.dir,'/Year',sep='')
# I commented out the above three input paths in the data.gen.r file since they are specified here and will change depending on scenarios

# new.dir <- "C:/Users/mddamian/Desktop/MARFIN MSE/Results/BSB std" # spare std folder for tricking loop past non convergent sims

dir.create(sce.dir)
dir.create(input.dir)
dir.create(input.season.dir)
dir.create(input.year.dir)

file.copy(from=paste(AM.code.dir,'/NSLSAP01.exe',sep=''),to=sce.dir) # copy the compiled .tpl file (executable file) to scenario specific folder
file.copy(from=paste(OM.code.dir,'/data.gen.r',sep=''), to=sce.dir) # two OM files added to scenario so it's easier for loop to call them
file.copy(from=paste(OM.code.dir,'/input.data.2021.r',sep=''), to=sce.dir) # reinitializes the starting conditions for each iteration below


n.iter = 200 # iterations (sampling variation)
n.AM = 10 # number of assessments
n_proj = 5 # number of years in an assessment cycle

# create matrices and arrays for each iteration to store outputs of interest for performance statistics
# most of these will need to be arrays to accomodate multiple scenarios
# save everything from report file for use later, e.g., adding performance metrics
pred.catch4d.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter))

obs.catch4d.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter))
pred.ssb.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
pred.survey.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
obs.survey.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
obs.ssb.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
pred.rec.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) 
obs.rec.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
pred.recdev.array <- array(NA, dim =c(n_t+n_proj*n.AM, n.AM+1, n.iter))
mean.rec <- array(NA, dim = c(n.AM+1, n.iter))
pred.survey.cpue.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
obs.survey.cpue.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter))
pred.abun.at.size.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_l, n.AM+1, n.iter))
obs.abun.at.size.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_l, n.AM+1, n.iter))
pred.comp.array <- array(NA, dim = c((n_t+n_proj*n.AM)*n_fleets, n_l, n.AM+1, n.iter))
obs.comp.array <- array(NA, dim = c((n_t+n_proj*n.AM)*n_fleets, n_l, n.AM+1, n.iter))
F.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter))
brp.array <- array(NA, dim = c(5, n.AM+1, n.iter))
bias.check <- array(NA, dim = c(n.iter))
noncon.counter = 0

# conv <- array(0, dim=c(n.AM+1, n.iter)) # assume 0 - converged, record as 1 below if it doesn't converge
iter = 1

### MSE LOOP ###

while(iter <= n.iter){
  
  source("input.data.2021.r") # reset initial conditions
  add.error = TRUE # whether to add sampling error - resetting initial conditions requires this to be set with each iteration
  source("data.gen.r") # run OM # error in weight_m dimensions 4/4/2022
  # run assessment model for the first time
  
  unlink("NSLSAP01.STD") # remove the previous .std file in the folder
  setwd(sce.dir)
  file.create("NSLSAP01.DAT")
  #makeAD("SSAP_M",verbose=F)
  runAD("NSLSAP01",verbose=T) # run AM

  fileNames=c("NSLSAP01.STD")
  if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
    # conv[1,iter] = 0 # record converged run
    file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
    # file.copy("NSLSAP01.STD", new.dir, overwrite = TRUE)
  }
  
  while(file.exists(fileNames)==FALSE){ # try this to see if we can avoid the spoof
    noncon.counter = noncon.counter + 1
    source("data.gen.R")
    runAD("NSLSAP01",verbose=T)
    if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
      # conv[1,iter] = 0 # record converged run
      file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
      # file.copy("NSLSAP01.STD", new.dir, overwrite = TRUE)
    }
  }
  # True popdy record 1
  obs.abun.at.size.array[1:n_t,,1,iter] <- PopDy$Abundance
  obs.ssb.array[1:n_t,1,iter] <- PopDy$SSB
  obs.rec.array[1:n_t,1,iter] <- PopDy$Recruits
  # else{
  #   if(file.exists(fileNames)==FALSE){
  #     conv[1,iter] = 1 #record non convergence
  #     # setwd(new.dir)
  #     # file.copy("NSLSAP01.STD", sce.dir, overwrite = TRUE) # "borrow" the last std file
  #     # setwd(sce.dir) # reset back to original working directory
  #     # file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy rep file to spoof "convergence"
  #     # do.AM = do.AM + 1 # and keep going
  #   }
  # }
  
  # figures for looking at model fits
  
  # set up projections
  do.AM = 0
  
  while(do.AM < n.AM){
    
    report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
    # get F reference points from the AM 
    fyear <- report$BeginYear
    lyear <- report$EndYear
    new_nt <- length(fyear:lyear)
    
    F_fleets <- t(report$Fishing_Mortality[,new_nt])
    F_sq <- sum(F_fleets) # F at status quo (current F), currently less than F0.1
    # F_msy <- report$FMSY # F at MSY (Scenario 1 of 2018 projections)
    # because we cannot use a functional SR relationship with BSB, we are treating F40%SPR as a proxy for Fmsy
    # F_30SPR <- report$F30SPR
    F_40spr <- report$F40SPR
    F_0.1 <- report$F0.1
    F_max <- report$Fmax
    # Implementation Error #1: random normal error added to Fmsy proxy
    # error based on mean and sd of apical F from assessment period
    f.err <- rnorm(1, mean = 0.05, sd = 0.05)
    # f.err
    F_msy <- F_40spr+f.err
    F_msy75 <- 0.75*F_msy # 75% of F at MSY, used for projections in BAM (Scenario 3 of 2018 projections)
    # F_P40 <- 0.94*F_msy # 94% of F at MSY to approximate F with P*=0.4 (Scenario 2 of 2018 projections)
    # F_P38 <- 0.92*F_msy # F with P* = 0.375, MP scoped by SSC in 2018, might have been what was actually used by the Council
    # F_msy60 <- 0.60*F_msy
    # F_msy50 <- 0.50*F_msy
    F_0 <- 0.001
    # F_max <- report$Fmax
    # F_0.1 <- report$F0.1 # F0.1
    Msy_scen <- c(F_msy,F_msy75,F_0,F_0.1,F_max) #updated, probably not worth testing those last two
    
    
    # using above F for projection 
    
    # Set up and run projection (PopDy - underlying population from OM)
    N_terminal = as.vector(tail(PopDy$Abundance,n=1))*exp(-(as.vector(tail(PopDy$F_3d,n=1))+as.vector(tail(PopDy$M_2d,n=1))))
    R_proj <- array(NA, dim = c(n_proj,1)) 
    R_devs_proj <- rnorm(n_proj, mean = 0, sd = sd_Rdevs) # based on last 10 years of assessment period
    R_devs_proj <- scale(R_devs_proj)
    Rec.mean <- mean(PopDy$Recruits[25:30])
    for(n in 1:n_proj){
      R_proj[n] <- Rec.mean*exp(R_devs_proj[n]) # don't forget to set to alpha if assuming avg rec conditions
    }
    F_apex_proj = rep(Msy_scen[5],n_proj)
    M_vec_proj = tail(M_m,n=1)
    G_proj = GM
    
    # If NOT using implementation error, use the following code:
    # tot.mean.F <- sum(mean(F_array[23:32,,1]),
    # mean(F_array[23:32,,2]),
    # mean(F_array[23:32,,3]))
    # 
    # mean(F_array[23:32,,1])/tot.mean.F # proportion from comm
    # mean(F_array[23:32,,2])/tot.mean.F # from rec
    # mean(F_array[23:32,,3])/tot.mean.F # from disc
    # F_fl_prop = c(0.08, 0.70, 0.21) # fixed average proportions of commercial, recreational and discard F based on avg of last 10 yrs
    
    # If using implementation error, use the following code:
    tot.mean.F <- sum(mean(F_array[23:32,,1]),
                      mean(F_array[23:32,,2]),
                      mean(F_array[23:32,,3]))
    F.assess <- F_array[23:32,,] # used later for implementation error #2
    F.assess <- F.assess/tot.mean.F
    F.assess.ten <- array(NA, dim = c(10,3))
    # logit transform proportions
    for(f in 1:3){
      for(i in 1:10){
        F.assess.ten[i,f] <- log(F.assess[i,f]/(1-F.assess[i,f]))
      }
    }
    F.assess.ten[3,2] <- log(1.432/abs(1-1.432)) # fix divide by - issue
    F.assess.ten/rowSums(F.assess.ten)
    # Draw two random variables for first and second fleet based on inverse logit-transforms
    x1 <- 1/(1+exp(-rnorm(1, mean = mean(F.assess.ten[,1], sd = sd(F.assess.ten[,1])))))
    p1 <- x1
    x2 <- 1/(1+exp(-rnorm(1, mean = mean(F.assess.ten[,2], sd = sd(F.assess.ten[,2])))))
    p2 <- x2*(1-p1)
    p3 <- 1-x1-x2*(1-x1)
    
    F_fl_prop<- c(p1, p2, p3)

    F_allocation_proj = array(NA,dim=c(n_proj,n_s,n_fleets)) # allocate F_apex_proj to fisheries - needs changing to reflect true allocation
    for(f in 1:n_fleets){
      F_allocation_proj[,,f] = F_fl_prop[f]*F_apex_proj
    }
    # to mimic the council management process for BSB, we need 3 scenarios: Fmsy, 0.75Fmsy and P*=0.4
    # To do the last one, we need to multiply the biomass variance by the value of P*
    # Important: we do not have a biomass sigma, so we will need to approximate the rate of Fmsy with P*=0.4
    # The SAFMC does not use linear HCRs that scale F with SSB, so status quo MPs are still constant catch scenarios
    # For the other MPs, adjusting selectivities makes the most sense for size-based measures
    # We can also plot estimated recreational and discard selectivities against effort to identify some optimal catch rate, 
    # e.g., the level of effort required to select for legal or larger fish
    
    # Need to think about where you run selectivity scenarios to approximate changes in size limits
    # Best way to do that is having scenarios made up of different selectivity parameter sets
    # Should probably go here before mse.proj is run, need to change what's used in mse.proj file and
    # maybe create a separate input.update file for selectivity runs
    
    #Finally, think about how to plot out selectivity vs. effort by rec and discard fleets to 
    #identify the necessary effort to select larger (legal) fish
    
    Stats.proj <- mse.proj(n_proj,N_terminal,R_proj,sels_4D,F_apex_proj,M_vec_proj,G_proj,F_allocation_proj)
    catch.proj.obj <- list(catch_4d = Stats.proj$Catch4d, catch_3d = Stats.proj$Catch3d)
    survey.proj.obj <- list(data.frame(Stats.proj$Survey))
    
    # Store original + projected values for performance statistics
    # set up arrays to store iteration averages by scenario
    # catch.matrix <- array(NA, dim = c(n.iter))
    # cv.matrix <- array(NA, dim = c(n.iter))
    # ssb.matrix <- array(NA, dim = c(n.iter))
    # 
    # catch.matrix <- 
    
    # update underlying population
    
    Abundance.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    Abundance.temp[,,1] = rbind(PopDy$Abundance[,,1],Stats.proj$N_proj[,,1])
    
    # Jie suggests calculating SSB here
    SSB.old <- as.matrix(PopDy$SSB)
    SSB.new <- as.matrix(Stats.proj$SSB_proj)
    SSB.temp <- rbind(SSB.old, SSB.new)
    
    for(i in 1:(n_t+n_proj)){
      
    }
    
    F.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    F.temp[,,1] = rbind(PopDy$F_3d[,,1],Stats.proj$F_proj) # calling F_3d from PopDy object in data.gen; stats.proj$F doesn't match up with input.data F_3d
    
    M.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    M.temp[,,1] = rbind(PopDy$M_2d[,,1],Stats.proj$M_proj)
    
    PopDy = list(Abundance = Abundance.temp, 
                 SSB=SSB.temp, # need more work to correctly compute
                 Recruits=c(PopDy$Recruits,Stats.proj$R_proj),F_3d=F.temp, M_2d=M.temp) # problem - F_3d changes after being affixed to PopDy
    
    # update model dimension and relevant quantities and generate new data sets (get ready for the subsequent assessment) 
    
    Year.proj.start = year_start + n_t
    Year.proj.end = Year.proj.start + n_proj -1 
    n_t = n_t + n_proj
    year_end_retro = Year.proj.end
    year_end = Year.proj.end
    
    source("C:/Users/matt.damiano/Desktop/MARFIN MSE/Operating Model and Projection code/input.data.2021.update.r")
    
    source("C:/Users/matt.damiano/Desktop/MARFIN MSE/Operating Model and Projection code/data.gen.proj.r")
    
    # run AM
    unlink("NSLSAP01.STD") # remove the previous .std file in the folder
    runAD("NSLSAP01",verbose=T) # run AM
    fileNames=c("NSLSAP01.STD")
    if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
      # conv[do.AM+1, iter] = 0
      file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
      do.AM = do.AM + 1
      # file.copy("NSLSAP01.STD", new.dir, overwrite = TRUE)
      # need to record convergence here, e.g. if STD file doesn't exist, conv = 1 (like in the popdy final)
      # but we also need to be able to restart from place if this occurs in the first assessment of an iteration - right now it will end the loop
      
    }
    while(file.exists(fileNames)==FALSE){ # try this to see if we can avoid the spoof
      noncon.counter = noncon.counter + 1
      ### Need to remove last three datapoints for nonconvergent runs
      # These include: Catch Data, Survey Data
      # survey data object is survey.data
      survey.data <- survey.data[1:(n_t-n_proj),] # removes last 3 datapoints if no convergence (seems to work 8-17-2022)
      # catch data object is catch.data, more complicated due to multiple fleets
      # needs to remove the last three data points for each fleet, could be done sequentially from the end
      catch.data <- catch.data[-((3*n_t-n_proj+1):(3*n_t)),] # removes last three from fleet three
      catch.data <- catch.data[-((2*n_t-n_proj+1):(2*n_t)),] # should remove last three of fleet two
      catch.data <- catch.data[-((n_t-n_proj+1):n_t),]
      source("C:/Users/matt.damiano/Desktop/MARFIN MSE/Operating Model and Projection code/data.gen.proj.r")
      runAD("NSLSAP01",verbose=T)
      if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
        # conv[1,iter] = 0 # record converged run
        file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
        do.AM = do.AM + 1
      }
    }
    # True popdy record 2
    obs.abun.at.size.array[1:(new_nt+n_proj),,do.AM+1,iter] <- PopDy$Abundance
    obs.ssb.array[1:(new_nt+n_proj),do.AM+1,iter] <- PopDy$SSB
    obs.rec.array[1:(new_nt+n_proj),do.AM+1,iter] <- PopDy$Recruits
  }
  # report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
  # fill matrices and arrays
  # catch3d.array[,iter] <- colSums(report$Catch_Pred) # Unknown error present: get an error with every other run that number of items not mult of replacement length
  # Store results from final report
  # pred.catch3d.array[1:(new_nt+n_proj),do.AM+1,iter] <- colSums(report$Catch_Pred)
  # conv.mat[iter,do.AM] <- con
  
  # Check bias
  # report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
  # report <- read.rep('Fmsy-2020-1.rep')
  # bias.check[iter] <- mean(apply(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100,2,mean))
  # bias.check[iter] <- max(abs(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100)) # maximum % bias of N-at-age
  iter = iter + 1 # update the iteration
} 

### LOOP TO RECORD MSE OUTPUTS ###

  iter = 1

  while(iter <= n.iter){ # must read in each report for each assessment over all iterations
    source("input.data.2021.r") # get starting n_t, year_start values
    do.AM = 0
    while(do.AM < n.AM){
      report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
      # get F reference points from the AM 
      fyear <- report$BeginYear
      lyear <- report$EndYear
      new_nt <- length(fyear:lyear)
      
      F_fleets <- t(report$Fishing_Mortality[,new_nt])
      F_sq <- sum(F_fleets) # F at status quo (current F), currently less than F0.1
      # F_msy <- report$FMSY # F at MSY (Scenario 1 of 2018 projections)
      # because we cannot use a functional SR relationship with BSB, we are treating F40%SPR as a proxy for Fmsy
      F_0.1 <- report$F0.1 # F0.1
      F_max <- report$Fmax
      F_40spr <- report$F40SPR
      F_msy <- F_spr40+f.err
      F_msy75 <- 0.75*F_msy # 75% of F at MSY, used for projections in BAM (Scenario 3 of 2018 projections)
      # F_P40 <- 0.94*F_msy # 94% of F at MSY to approximate F with P*=0.4 (Scenario 2 of 2018 projections)
      # F_P38 <- 0.92*F_msy # F with P* = 0.375, MP scoped by SSC in 2018, might have been what was actually used by the Council
      # F_max <- report$Fmax
      Msy_scen <- Msy_scen <- c(F_msy,F_msy75,F_0,F_0.1,F_max) #updated, probably not worth testing those last two
      
      
      # using above F for projection 
      
      # Set up and run projection (PopDy - underlying population from OM)
      # N_terminal = as.vector(tail(PopDy$Abundance,n=1))*exp(-(as.vector(tail(PopDy$F_3d,n=1))+as.vector(tail(PopDy$M_2d,n=1))))
      # R_proj = rep(alpha,n_proj) # need to be changed later, currently set at mean recruitment without variation from year to year
      # F_apex_proj = rep(Msy_scen[4],n_proj)
      # M_vec_proj = tail(M_m,n=1)
      # G_proj = GM
      # F_fl_prop = c(0.17, 0.31, 0.52) # fixed average proportions of commercial, recreational and discard F based on avg of last 10 yrs
      # # updated 6-23-22, original proportions were wrong. 
      # F_allocation_proj = array(NA,dim=c(n_proj,n_s,n_fleets)) # allocate F_apex_proj to fisheries - needs changing to reflect true allocation
      # for(f in 1:n_fleets){
      #   F_allocation_proj[,,f] = F_fl_prop[f]*F_apex_proj
      # }
      # 
      Year.proj.start = year_start + n_t
      Year.proj.end = Year.proj.start + n_proj -1 
      n_t = n_t + n_proj
      year_end_retro = Year.proj.end
      year_end = Year.proj.end
      
      # record results from first assessment (fit to OM)
      pred.catch4d.array[1:new_nt,1,do.AM+1,iter] <- report$Catch_Pred[1,]
      pred.catch4d.array[1:new_nt,2,do.AM+1,iter] <- report$Catch_Pred[2,]
      pred.catch4d.array[1:new_nt,3,do.AM+1,iter] <- report$Catch_Pred[3,]
      obs.catch4d.array[1:new_nt,1,do.AM+1,iter] <- report$Catch_Obs[1,]
      obs.catch4d.array[1:new_nt,2,do.AM+1,iter] <- report$Catch_Obs[2,]
      obs.catch4d.array[1:new_nt,3,do.AM+1,iter] <- report$Catch_Obs[3,]
      # survey
      pred.survey.array[1:new_nt,do.AM+1,iter] <- report$Survey_Index_Pred
      # obs.survey.array[1:new_nt,do.AM+1,iter] <- report$Survey_Index_Obs
      # ssb and recruitment
      # pred.ssb.array[1:new_nt, do.AM+1, iter] <- report$Spawning_stock_Biomass
      pred.ssb.array[1:new_nt, do.AM+1, iter] <- report$Spawning_stock_Biomass_input
      pred.rec.array[1:new_nt, do.AM+1, iter] <- report$recruitment_Pred
      mean.rec[do.AM+1, iter] <- report$`Mean_Recruitment:`
      pred.recdev.array[1:new_nt,do.AM+1, iter] <- report$recruitment_log_Dev
      # numbers at size and catch composition
      pred.abun.at.size.array[1:new_nt, , do.AM+1, iter] <- report$Abundance_at_Size
      pred.comp.array[1:(new_nt*n_fleets), , do.AM+1, iter] <- report$Catch_Comp_Pred
      obs.comp.array[1:(new_nt*n_fleets), , do.AM+1, iter] <- report$Catch_Comp_Obs
      # F
      F.array[1:new_nt, 1, do.AM+1, iter] <- report$Fishing_Mortality[1,]
      F.array[1:new_nt, 2, do.AM+1, iter] <- report$Fishing_Mortality[2,]
      F.array[1:new_nt, 3, do.AM+1, iter] <- report$Fishing_Mortality[3,]
      brp.array[1:5, do.AM+1, iter] <- Msy_scen
      source("C:/Users/matt.damiano/Desktop/MARFIN MSE/Operating Model and Projection code/input.data.2021.update.r")
      do.AM = do.AM+1
    }
    # call the most recent report file
    report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
    # fill matrices and arrays
    # catch3d.array[,iter] <- colSums(report$Catch_Pred) # Unknown error present: get an error with every other run that number of items not mult of replacement length
    # report <- read.rep(paste(scenario,'-',2017,'-',1,".","rep",sep=""))
    
    pred.survey.array[1:(new_nt+n_proj),do.AM+1,iter] <- report$Survey_Index_Pred
    obs.survey.array[1:(new_nt+n_proj),do.AM+1,iter] <- report$Survey_Index_Obs
    # catch 4d
    pred.catch4d.array[1:(new_nt+n_proj),1,do.AM+1,iter] <- report$Catch_Pred[1,]
    pred.catch4d.array[1:(new_nt+n_proj),2,do.AM+1,iter] <- report$Catch_Pred[2,]
    pred.catch4d.array[1:(new_nt+n_proj),3,do.AM+1,iter] <- report$Catch_Pred[3,]
    obs.catch4d.array[1:(new_nt+n_proj),1,do.AM+1,iter] <- report$Catch_Obs[1,]
    obs.catch4d.array[1:(new_nt+n_proj),2,do.AM+1,iter] <- report$Catch_Obs[2,]
    obs.catch4d.array[1:(new_nt+n_proj),3,do.AM+1,iter] <- report$Catch_Obs[3,]
    # ssb and recruitment
    # pred.ssb.array[1:(new_nt+n_proj), do.AM+1, iter] <- report$Spawning_stock_Biomass
    pred.ssb.array[1:(new_nt+n_proj), do.AM+1, iter] <- report$Spawning_stock_Biomass_input
    pred.rec.array[1:(new_nt+n_proj), do.AM+1, iter] <- report$recruitment_Pred
    mean.rec[do.AM+1, iter] <- report$`Mean_Recruitment:`
    pred.recdev.array[1:(new_nt+n_proj),do.AM+1, iter] <- report$recruitment_log_Dev
    # numbers at size and catch composition
    pred.abun.at.size.array[1:(new_nt+n_proj), , do.AM+1, iter] <- report$Abundance_at_Size
    pred.comp.array[1:((new_nt+n_proj)*n_fleets), , do.AM+1, iter] <- report$Catch_Comp_Pred
    obs.comp.array[1:((new_nt+n_proj)*n_fleets), , do.AM+1, iter] <- report$Catch_Comp_Obs
    # F
    F.array[1:(new_nt+n_proj), 1, do.AM+1, iter] <- report$Fishing_Mortality[1,]
    F.array[1:(new_nt+n_proj), 2, do.AM+1, iter] <- report$Fishing_Mortality[2,]
    F.array[1:(new_nt+n_proj), 3, do.AM+1, iter] <- report$Fishing_Mortality[3,]
    brp.array[1:5, do.AM+1, iter] <- Msy_scen
    
    # Check bias
    # report <- read.rep('Fmsy-2020-1.rep')
    bias.check[iter] <- mean(apply(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100,2,mean))
    # bias.check[iter] <- max(abs(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100)) # maximum % bias of N-at-age
    iter = iter + 1 # update the iteration
  }
  
  # calculate average annual variability for commercial fleet
  AAV <- array(NA, dim=c(n_fleets,n_t-1,n.AM+1,n.iter))
  for(i in 1:n.iter){
    for(a in 1:(n.AM+1)){
      for(t in 1:(n_t-1)){
        for(f in 1:n_fleets){
          AAV[f,t,a,i] <- pred.catch4d.array[t+1,f,a,i] - pred.catch4d.array[t,f,a,i]
        }
      }
    }
  }
  
  # We want the last year because we are measuring performance of MPs over long-term;
  # presenting uncertainty in prior assessment cycles doesn't tell us what we want
  # Final assessment results
  pred.catch.comm.final <- pred.catch4d.array[,1,(n.AM+1),]
  pred.catch.rec.final <- pred.catch4d.array[,2,(n.AM+1),]
  pred.catch.disc.final <- pred.catch4d.array[,3,(n.AM+1),]
  obs.catch.comm.final <- obs.catch4d.array[,1,(n.AM+1),]
  obs.catch.rec.final <- obs.catch4d.array[,2,(n.AM+1),]
  obs.catch.disc.final <- obs.catch4d.array[,3,(n.AM+1),]
  pred.ssb.final <- pred.ssb.array[,(n.AM+1),]
  obs.ssb.final <- obs.ssb.array[,(n.AM+1),]
  pred.rec.final <- pred.rec.array[,(n.AM+1),]
  pred.recdev.final <- pred.recdev.array[,(n.AM+1),]
  obs.rec.final <- obs.rec.array[,(n.AM+1),]
  pred.rec.first <- pred.rec.array[,1,]
  obs.rec.first <- obs.rec.array[,1,]
  abd.size.final <- pred.abun.at.size.array[,,(n.AM+1),]
  pred.survey.final <- pred.survey.array[,(n.AM+1),]
  obs.survey.final <- obs.survey.array[,(n.AM+1),]
  
  
  # We also want exploitation rates so we can calculate some performance measures consistent with Bohaboy et al. 2022
  # need N summed over size classes 
  abd.3d <- apply(abd.size.final, c(1,3), sum)
  exp.comm.final <- pred.catch.comm.final/abd.3d # commercial exploitation rates
  exp.rec.final <- pred.catch.rec.final/abd.3d # recreational exploitation rates
  exp.disc.final <- pred.catch.disc.final/abd.3d
  # catch rate proxy: rec catch/rec exp rate
  catch.rate <- pred.catch.rec.final/exp.rec.final
  # number of legal sized fish in the population
  abd.legal <- abd.size.final[,8:17,] # legal size classes for 11 inch size limit
  abd.legal.3d <- apply(abd.legal, c(1,3), sum)
  abd.legal.prop <- abd.legal.3d/abd.3d
  
  true.abd.size.final <- obs.abun.at.size.array[,,(n.AM+1),]
  true.abd.3d <- apply(true.abd.size.final, c(1,3), sum)
  
  # Convergence rate
  conv <- ((n.AM+1)*n.iter)/(((n.AM+1)*n.iter)+noncon.counter)
  
  # pred.catch3d.means <- rowMeans(pred.catch3d.array[,1:(n.AM+1),1:n.iter], dims = 2, na.rm=T)
  # pred.catch3d.means <- apply(pred.catch3d.array,c(1,3),mean,na.rm=T)
  # # obs.catch3d.means <- rowMeans(obs.catch3d.array[,1:(n.AM+1),1:n.iter], dims = 2, na.rm=T)
  # obs.catch3d.means <- apply(obs.catch3d.array,c(1,3),mean,na.rm=T)
  
  # Store results to save and load later as needed
  save(scenario,obs.catch4d.array,pred.catch4d.array,obs.survey.array,pred.survey.array,obs.abun.at.size.array,pred.abun.at.size.array,
       pred.ssb.array,obs.ssb.array,pred.rec.array,mean.rec,pred.comp.array,obs.comp.array,F.array,brp.array,pred.catch.comm.final,
       pred.catch.rec.final,pred.catch.disc.final,obs.catch.comm.final,obs.catch.rec.final,obs.catch.disc.final,pred.ssb.final,obs.ssb.final,
       pred.rec.final,obs.rec.final,abd.size.final,pred.survey.final,obs.survey.final,abd.3d,exp.comm.final,exp.rec.final,exp.disc.final,
       catch.rate,abd.legal,abd.legal.3d,abd.legal.prop,true.abd.size.final, true.abd.3d,pred.recdev.array,pred.recdev.final,
       AAV, bias.check, conv,
       file = "Results.Rdata", envir = .GlobalEnv)
  
#### END MSE #####
  

#### VISUALS ####
  
# load file as needed for visuals
  rm(list=ls())
  setwd('C:/Users/matt.damiano/Desktop/MARFIN MSE/BSB Results')
  # load('C:/Users/mddamian/Desktop/MARFIN MSE/Results/Final Sims/BSB P38 mean rec est sd = 100 8-24-2022.Rdata')
  # report <- 
  load('BSB FMSY SL rec rec.Rdata')
  # Diagnostic plots
  year_start = 1990
  year_end = 2071
  n.iter = 200
  n.AM = 10
  n_t = 82
  size_bm = c(90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490)
  
  # true.abd.size.final <- obs.abun.at.size.array[,,(n.AM+1),]
  # true.abd.3d <- apply(true.abd.size.final, c(1,3), sum)

  # final proportion of abundance
  # equil. abd = approx 8e7
  # eq.abd.size <- pred.abun.at.size.array[78,,18,50]
  # final.abd.prop <- (eq.abd.size/(8e+07))*100
  # final.abd.prop
  # # Jitter plots for time series (SSB, survey CPUE, total catch)
  # # 7-14-2022: not sure how to jitter over both AM and iterations...using just the averages for now
  # 
  # # SSB and Rec, obs and pred
  # # averaging done here are over assessments, so we end up with one mean that has been averaged over iterations first, then assessments
  # # Try averaging over assessments first, then you can actually jitter over iterations - this will help us see where the variation is occurring: assessment? iteration?
  # par(mfrow = c(1,3))
  # par(mfrow = c(1,1))
  # plot(seq(year_start,year_end,1),rowMeans(obs.ssb.final, dim = 1), type="l", ylim=c(0,3.5e10), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="True SSB (eggs)", xlab="Year",main=scenario)
  # for (i in 1:n.iter){
  #   points(seq(year_start,year_end,1),jitter(obs.ssb.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.ssb.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  # }
  # 
  # Pred Recruitment
  
  ssb.means <- rowMeans(pred.ssb.final, dim =1)
  MSST = ssb.means[27]/1.15
  MSST = array(MSST, dim = c(79,1))
  
  # need to fix an error
  obs.ssb.final[1,1:200] <- 13899880922

  # True recruitment (OM from first two size bins)
  # plot(seq(year_start,year_end,1),rowMeans(obs.rec.final, dim=1), type="l",ylim=c(0,8e7), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla=" True Recruits (to size bins 1 and 2)", xlab="Year")
  # for (i in 1:n.iter){
  #   points(seq(year_start,year_end,1),jitter(obs.rec.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.rec.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  # }
  # Recruitment
  windows(width=15,height=10)
  par(mfrow = c(1,3))
  plot(seq(year_start,year_end,1),rowMeans(pred.rec.final, dim=1), type="l",ylim=c(0,max(pred.rec.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Recruits", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.rec.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.rec.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
    #points(seq(year_start,year_end,1),jitter(obs.rec.final[,i], factor=4), type="l", col= "black", lwd=2)
  }
  # Abundance
  plot(seq(year_start,year_end,1),rowMeans(abd.3d, dim=1), type="l",ylim=c(0,max(abd.3d)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Abundance", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(abd.3d[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(abd.3d[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
    #points(seq(year_start,year_end,1),jitter(true.abd.3d[,i], factor=4), type="l", col= "black", lwd=2)
  }
  # SSB
  plot(seq(year_start,year_end,1),rowMeans(pred.ssb.final, dim = 1), type="l", ylim=c(0,max(pred.ssb.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="SSB (eggs)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.ssb.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.ssb.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
    #points(seq(year_start,year_end,1),jitter(obs.ssb.final[,i], factor=4), type="l", col= "black", lwd=2)
    # points(seq(year_start,year_end,1),MSST, type="p", col= "black", lwd=0.2)
  }
  # # True abundance 
  # plot(seq(year_start,year_end,1),rowMeans(true.abd.3d, dim=1), type="l",ylim=c(0,1.6e8), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla=" True total N", xlab="Year")
  # for (i in 1:n.iter){
  #   points(seq(year_start,year_end,1),jitter(true.abd.3d[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(true.abd.3d[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  # }
  
  
  # Look at catch by fleet
  windows(width=15,height=10)
  par(mfrow = c(3,2))
  plot(seq(year_start,year_end,1),rowMeans(pred.catch.rec.final, dim =1), type="l",ylim=c(0,max(pred.catch.rec.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Rec Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.catch.rec.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.catch.rec.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(obs.catch.rec.final, dim =1), type="l",ylim=c(0,max(pred.catch.rec.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="True Rec Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(obs.catch.rec.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.catch.rec.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(pred.catch.comm.final, dim =1), type="l",ylim=c(0,max(pred.catch.comm.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Comm Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.catch.comm.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.catch.comm.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(obs.catch.comm.final, dim =1), type="l",ylim=c(0,max(pred.catch.comm.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="True Comm Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(obs.catch.comm.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.catch.comm.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(pred.catch.disc.final, dim =1), type="l",ylim=c(0,max(pred.catch.disc.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Disc Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.catch.disc.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.catch.disc.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(obs.catch.disc.final, dim =1), type="l",ylim=c(0,max(pred.catch.disc.final)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="True Disc Catch (fish)", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(obs.catch.disc.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.catch.disc.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  
  windows(width=15,height=10)
  par(mfrow = c(1,2))
  # Survey CPUE Pred and Obs
  plot(seq(year_start,year_end,1),rowMeans(pred.survey.final, dim =1), type="l",ylim=c(0,200), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Survey CPUE", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(pred.survey.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.survey.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  plot(seq(year_start,year_end,1),rowMeans(obs.survey.final, dim=1), type="l",ylim=c(0,200), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="True Survey CPUE", xlab="Year")
  for (i in 1:n.iter){
    points(seq(year_start,year_end,1),jitter(obs.survey.final[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.survey.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  }
  # Total catch
  # plot(seq(year_start,year_end,1),rowMeans(pred.catch3d.means, dim=1), type="l",ylim=c(0,1e+06), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Catch (fish)", xlab="Year")
  # for (i in 1:n.iter){
  #   points(seq(year_start,year_end,1),jitter(pred.catch3d.means[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.catch3d.means[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  # }
  # plot(seq(year_start,year_end,1),rowMeans(obs.catch3d.means, dim=1), type="l",ylim=c(0,1e+06), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Observed Catch (fish)", xlab="Year")
  # for (i in 1:n.iter){
  #   points(seq(year_start,year_end,1),jitter(obs.catch3d.means[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(obs.catch3d.means[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
  # }
  # We also want box and whisker plots for terminal year estimates to see uncertainty there
  
  # old.par <- par(mar = c(0,0,0,0))
  # par(old.par)
  
  # Abundance at size bubble plot
  abd.size.final <- abd.size.final[1:n_t,,] # need to remove final two years so uncertain recruitment doesn't mute size structure
  abd.size.final.mean <- rowMeans(abd.size.final, dim = 2) # takes the mean abundance at size over 200 iterations of the final assessment
  
  # abd.size.means <- rowMeans(pred.abun.at.size.array[,,1:(n.AM+1),1:n.iter], dim=2, na.rm=T)
  # abd.size.means <- apply(pred.abun.at.size.array,c(1,2),mean,na.rm=T) # take mean over AMs + iterations
  # Create new matrix of averages over iterations to use plot function with
  # AvgAbdSize <- matrix(NA, nrow = n_t, ncol = n_l)
  # AvgAbdSize <- rowMeans(abun.at.size.array, dims = 2) # average numbers at length over all sims
  # 
  # Run this or it won't work:
  #Bubble function
  #######################################################################
  panel.bubble.for.R<-function(x,y,z,subscripts,scalem=4,...){
    #panel.grid(col="lightgrey",lwd=1,h=-1,v=-1)
    cex<-abs(z) # size of points depends on the absolute value of z and is consistent across panels
    cex<-scalem*sqrt(cex/max(cex))
    flagplus<- z>=0
    panel.xyplot(x[flagplus],y[flagplus],pch=16,col="gray",cex=cex[flagplus])
    panel.xyplot(x[!flagplus],y[!flagplus],pch=1,col="gray",cex=cex[!flagplus])
  }
  
  
  PlotAbun=function()
  {
    library(lattice)
    Abun=abd.size.final.mean
    Season=1 # no seasonality in model
    Year=seq(year_start,year_end,1)
    sizebins=size_bm # had to change to midpoints to match abundance matrix
    
    AbunB=matrix(NA,ncol=3,nrow=length(Year)*length(sizebins))
    itemp=1
    
    for (i in Year)
    {
      iL=i-year_start+1
      
      for (j in 1:length(sizebins))
      {
        if (Season>1)
        {
          AbunB[itemp,]=c(i,sizebins[j],Abun[(iL-1)*Season+1,j])
          itemp=itemp+1
        }else
        {
          AbunB[itemp,]=c(i,sizebins[j],Abun[iL,j])
          itemp=itemp+1
        }
        
      }
    }
    AbunB<-as.data.frame(AbunB)    
    names(AbunB)=c("Year","CL","Abun")
    
    plot<-xyplot(CL ~ Year,data=AbunB,layout=c(1,1),par.strip.text=list(cex=1.3),
                 as.table=T,subscript=T,z=AbunB[,3],
                 # main=paste(scenario),
                 scales=list(y=list(cex=1.3),x=list(cex=1.3)),panel=panel.bubble.for.R,scalem=4,
                 ylab=list(label="Numbers at Length (mm)",cex=1.3),xlab=list(label="Year", cex=1.3))
    print(plot)  
    #  par(mfrow=c(1,1)
  }
  windows(width=15,height=10)
  PlotAbun()
  
  
  # Histrograms of mean catch by fleet
  windows(width=15,height=10)
  par(mfrow = c(1,3))
  # commercial
  hist(apply(pred.catch4d.array[1:n_t, 1, 1:(n.AM+1),1:n.iter],1,mean, na.rm=T),xlab='Mean Catch (commercial)',main="Comm")
  abline(v = median(apply(pred.catch4d.array[1:n_t,1,1:(n.AM+1),1:n.iter],1,mean), na.rm=T), col = "red",lwd=3)
  # recreational
  hist(apply(pred.catch4d.array[1:n_t, 2, 1:(n.AM+1),1:n.iter],1,mean, na.rm=T),xlab='Mean Catch (recreational)',main="Rec")
  abline(v = median(apply(pred.catch4d.array[1:n_t,2,1:(n.AM+1),1:n.iter],1,mean), na.rm=T), col = "red",lwd=3)
  # discard
  hist(apply(pred.catch4d.array[1:n_t, 3, 1:(n.AM+1),1:n.iter],1,mean, na.rm=T),xlab='Mean Catch (discard)',main="Disc")
  abline(v = median(apply(pred.catch4d.array[1:n_t,3,1:(n.AM+1),1:n.iter],1,mean), na.rm=T), col = "red",lwd=3)
  # hist(bias.check[1:n.iter], xlab = 'Percent bias', main = "Bias check")
  
  # # Only needs to be re-run if size limit is lowered
  # abd.legal <- abd.size.final[,8:17,] # legal size classes for 11 inch size limit
  # abd.legal.3d <- apply(abd.legal, c(1,3), sum)
  # abd.legal.prop <- abd.legal.3d/abd.3d
  
  # Exploitation rate (season length) and catch rate
  windows(width=15,height=10)
  par(mfrow = c(1,3))
  hist(apply(exp.rec.final[1:n_t,1:n.iter],1,mean), xlab = "Mean exploitation rate (recreational)", main = "Exp rate/Season length")
  abline(v = median(apply(exp.rec.final[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)
  hist(apply(catch.rate[1:n_t,1:n.iter],1,mean), xlab = "Mean catch rate (recreational)", main = "Catch rate")
  abline(v = median(apply(catch.rate[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)
  hist(apply(abd.legal.prop[1:n_t,1:n.iter],1,mean), xlab = "Mean proportion of legal-sized fish (recreational)", main = "Proportion legal")
  abline(v = median(apply(abd.legal.prop[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)
  
  # reset graphing parameters
  # par(mfrow=c(1,2))
  # # plot average annual variability for commercial fleet
  # AAV.means <- array(NA, dim = c(n_t-1,n.AM+1,n.iter))
  # AAV.means <- rowMeans(AAV[1,,1:(n.AM+1),1:n.iter], dim=2, na.rm=T)
  # # seqAvgAAV <- rowMeans(AAV, dims = 2) # this is the actual AAV
  # # Year=seq(report$BeginYear,report$EndYear,1)
  # plot(seq(year_start,year_end-1,1),rowMeans(AAV.means[,1:(n.AM+1)], na.rm=T), type = "l", xlab = "Year", ylab = "AAV (number of fish)", main = "Average Annual Variability (Commercial)",
  #      lwd = 3)
  # 

# Recruitment from first assessment diagnostics (if needed)
#   windows(width=15,height=10)
#   par(mfrow = c(1,1))
#   plot(seq(year_start,year_end,1),rowMeans(pred.rec.first, na.rm=T), type="l",ylim=c(0,max(pred.rec.array[,1,], na.rm=T)), col=rgb(160,32,240,255,maxColorValue=255), lwd=3, yla="Predicted Recruits (to size bins 1 and 2)", xlab="Year")
#     points(seq(year_start,year_end,1),jitter(pred.rec.first[,i], factor=4), type="l", col=rgb(160,32,240,(255-abs(round(dim(pred.rec.final[,1:n.iter])[2]/2)-i))/9,maxColorValue=255), lwd=2)
#     points(seq(year_start,year_end,1),rowMeans(obs.rec.first, na.rm=T), type="l", col= "black", lwd=2)
# par(mfrow=c(1,1))
# plot(pred.recdev.array[,18,100], main = "Predicted rec devs with OM recdevs centered")
# points(pred.recdev.array[,18,100], col = "red")
# legend("topright", legend = c("iter 1", "iter 2"),
#        col=c("black", "red"), lty=3)

# boxplots
# need to transpose these for them to give us boxplots by year
t.pred.rec.final <- t(pred.rec.final)
t.obs.rec.final <- t(obs.rec.final)
t.pred.rec.first <- t(pred.rec.first)
obs.rec.first <- obs.rec.array[,1,]
t.obs.rec.first <- t(obs.rec.first)
# 
# boxplot(pred.rec.final[79,]) #bias in final year
# boxplot(t.pred.rec.final) # final assessment
# boxplot(t.pred.rec.first) # first assessment
windows(width=15,height=10)
boxplot((t.pred.rec.final-t.obs.rec.final)/t.obs.rec.final, ylim=c(-1,1)) # deviation from true recruitment

boxplot((t.pred.rec.first-t.obs.rec.first)/t.obs.rec.first, ylim=c(-1,1)) # deviation from true recruitment

# Quick checks for tradeoff reference points
median(apply(pred.catch4d.array[1:n_t,1,1:(n.AM+1),1:n.iter],1,mean), na.rm=T)
median(apply(pred.catch4d.array[1:n_t,2,1:(n.AM+1),1:n.iter],1,mean), na.rm=T)
median(apply(exp.rec.final[1:n_t, 1:n.iter],1,mean))
median(apply(catch.rate[1:n_t, 1:n.iter],1,mean))
median(apply(abd.legal.prop[1:n_t, 1:n.iter],1,mean))

# P38 reference
# median(apply(pred.catch4d.array[1:n_t,1,1:(n.AM+1),1:n.iter],1,mean), na.rm=T)
# [1] 1833354
# > median(apply(pred.catch4d.array[1:n_t,2,1:(n.AM+1),1:n.iter],1,mean), na.rm=T)
# [1] 2462357
# > median(apply(exp.rec.final[1:n_t, 1:n.iter],1,mean))
# [1] 0.03254346
# > median(apply(catch.rate[1:n_t, 1:n.iter],1,mean))
# [1] 81592915
# > median(apply(abd.legal.prop[1:n_t, 1:n.iter],1,mean))
# [1] 0.0528663
# > 
get.perdif <- function(x,y){
  perdif <- ((x-y)/y)*100
  return(perdif)
}
get.perdif(0.1232371,0.0528663)




