# Size-structured MSE 
# MSE projection file
# M. Damiano
# Last updated: 8-23-2023

# load libraries for running AD Model Builder code
require(PBSadmb)
require(psychTools) 

# Set up your Directories for the assessment model, the operating model, and results
# You will need to set up all directories before you can start using this model
AM.code.dir = "C:..."
OM.code.dir = "C:..."
result.dir = "C:..."

# Remember to read in the reptoRlist.R file to connect with ADMB
source('C:...reptoRlist.r')

# condition the population dynamics - OM set-up
setwd(OM.code.dir)
source("Functions.r")
source("input.data.r")
source("set.pars.r")
source("mse.proj.r")

########################################## simulation paths setting #######################################
scenario         = c('') # give scenario a unique name

# Set directories
sce.dir          = paste(result.dir,'/',scenario,sep='')
input.dir        = paste(result.dir,'/',scenario,'/InputFiles',sep='')
input.season.dir = paste(input.dir,'/Season',sep='')
input.year.dir   = paste(input.dir,'/Year',sep='')

dir.create(sce.dir)
dir.create(input.dir)
dir.create(input.season.dir)
dir.create(input.year.dir)

file.copy(from=paste(AM.code.dir,'/NSLSAP01.exe',sep=''),to=sce.dir) # copy the compiled .tpl file (executable file) to scenario specific folder
file.copy(from=paste(OM.code.dir,'/data.gen.r',sep=''), to=sce.dir) # two OM files added to scenario so it's easier for loop to call them
file.copy(from=paste(OM.code.dir,'/input.data.r',sep=''), to=sce.dir) # reinitializes the starting conditions for each iteration below


n.iter = 200 # iterations (sampling variation)
n.AM = 10 # number of assessments
n_proj = 5 # number of years in an assessment cycle

# create matrices and arrays for each iteration to store outputs of interest for performance statistics
# save everything from report file for use later, e.g., adding performance metrics
# Note: "observed" or "obs." arrays store values from the operating model: our "truth"
pred.catch4d.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter)) # predicted catch by year, AM, and iteration
obs.catch4d.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter)) # observed catch
pred.ssb.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # predicted ssb
pred.survey.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # predicted survey 
obs.survey.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # observed survey
obs.ssb.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # observed ssb
pred.rec.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # predicted recruitment
obs.rec.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # observed recruitment
pred.recdev.array <- array(NA, dim =c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # predicted recruitment deviations
mean.rec <- array(NA, dim = c(n.AM+1, n.iter)) # total mean recruitment
pred.survey.cpue.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # predicted survey CPUE
obs.survey.cpue.array <- array(NA, dim = c(n_t+n_proj*n.AM, n.AM+1, n.iter)) # observed survey CPUE
pred.abun.at.size.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_l, n.AM+1, n.iter)) # predicted N at size
obs.abun.at.size.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_l, n.AM+1, n.iter)) # observed N at size
pred.comp.array <- array(NA, dim = c((n_t+n_proj*n.AM)*n_fleets, n_l, n.AM+1, n.iter)) # predicted length comps (catch)
obs.comp.array <- array(NA, dim = c((n_t+n_proj*n.AM)*n_fleets, n_l, n.AM+1, n.iter)) # observed length comps
F.array <- array(NA, dim = c(n_t+n_proj*n.AM, n_fleets, n.AM+1, n.iter)) # F by fleet for each AM and iteration
brp.array <- array(NA, dim = c(5, n.AM+1, n.iter)) # stores biological reference points, e.g., FMSY
bias.check <- array(NA, dim = c(n.iter)) # for checking bias in fit
noncon.counter = 0 # for counting non-convergent runs

# starting iteration
iter = 1

### MSE LOOP ###

while(iter <= n.iter){
  
  source("input.data.r") # reset initial conditions
  add.error = TRUE # whether to add sampling error - resetting initial conditions requires this to be set with each iteration
  source("data.gen.r") # run OM 
  # run assessment model for the first time
  unlink("NSLSAP01.STD") # remove the previous .std file in the folder
  setwd(sce.dir)
  file.create("NSLSAP01.DAT")
  runAD("NSLSAP01",verbose=T) # run AM

  fileNames=c("NSLSAP01.STD")
  if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
    file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
  }
  
  # These next few lines of code are to record any non-convergence, kick the run out and start over
  while(file.exists(fileNames)==FALSE){ 
    noncon.counter = noncon.counter + 1
    source("data.gen.R")
    runAD("NSLSAP01",verbose=T)
    if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
      file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
    }
  }
  # Record "true" population dynamics from the OM
  obs.abun.at.size.array[1:n_t,,1,iter] <- PopDy$Abundance
  obs.ssb.array[1:n_t,1,iter] <- PopDy$SSB
  obs.rec.array[1:n_t,1,iter] <- PopDy$Recruits

  # set up projections
  do.AM = 0
  
  while(do.AM < n.AM){
    
    report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
    # get F reference points from the AM 
    fyear <- report$BeginYear
    lyear <- report$EndYear
    new_nt <- length(fyear:lyear)
    F_fleets <- t(report$Fishing_Mortality[,new_nt])
    F_40spr <- report$F40SPR
    F_0.1 <- report$F0.1
    F_max <- report$Fmax
    # Implementation Error #1: random normal error added to Fmsy proxy
    # error based on mean and sd of apical F from assessment period
    f.err <- rnorm(1, mean = 0.05, sd = 0.05)
    F_msy <- F_40spr+f.err
    F_msy75 <- 0.75*F_msy # 75% of F at MSY, used in SEDAR projections
    F_0 <- 0.001
    Msy_scen <- c(F_msy,F_msy75,F_0,F_0.1,F_max) # Put whatever scenarios by BRP you want to test in this vector
    
    # using above F for projection 
    
    # Set up and run projection (PopDy - underlying population from OM)
    N_terminal = as.vector(tail(PopDy$Abundance,n=1))*exp(-(as.vector(tail(PopDy$F_3d,n=1))+as.vector(tail(PopDy$M_2d,n=1))))
    R_proj <- array(NA, dim = c(n_proj,1)) 
    R_devs_proj <- rnorm(n_proj, mean = 0, sd = sd_Rdevs) # recdevs based on last 10 years of assessment period
    R_devs_proj <- scale(R_devs_proj)
    Rec.mean <- mean(PopDy$Recruits[25:30]) # this bases the projected mean recruitment based on a recent avg - you will want to change this AND the dimensions
    for(n in 1:n_proj){
      R_proj[n] <- Rec.mean*exp(R_devs_proj[n]) # don't forget to set to alpha if assuming avg rec conditions
    }
    F_apex_proj = rep(Msy_scen[5],n_proj) # This is where you select your BRP scenario, e.g., Msy_scen[5] is whatever is 5th in that vector
    M_vec_proj = tail(M_m,n=1)
    G_proj = GM
    
    # If NOT using allocation implementation error, use the following code:
    # tot.mean.F <- sum(mean(F_array[23:32,,1]),
    # mean(F_array[23:32,,2]),
    # mean(F_array[23:32,,3]))
    # 
    # mean(F_array[23:32,,1])/tot.mean.F # proportion from comm
    # mean(F_array[23:32,,2])/tot.mean.F # from rec
    # mean(F_array[23:32,,3])/tot.mean.F # from disc
    # F_fl_prop = c(0.08, 0.70, 0.21) # fixed average proportions of commercial, recreational and discard F based on avg of last 10 yrs
    
    # If using allocation implementation error, use the following code:
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
  
    # update underlying population
    
    Abundance.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    Abundance.temp[,,1] = rbind(PopDy$Abundance[,,1],Stats.proj$N_proj[,,1])
    
    # update SSB
    SSB.old <- as.matrix(PopDy$SSB)
    SSB.new <- as.matrix(Stats.proj$SSB_proj)
    SSB.temp <- rbind(SSB.old, SSB.new)
    
    # Update F
    F.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    F.temp[,,1] = rbind(PopDy$F_3d[,,1],Stats.proj$F_proj) # calling F_3d from PopDy object in data.gen; stats.proj$F doesn't match up with input.data F_3d
    
    # Update M
    M.temp = array(NA, c(n_proj+n_t,n_l,n_s))
    M.temp[,,1] = rbind(PopDy$M_2d[,,1],Stats.proj$M_proj)
    
    # Put all in one list
    PopDy = list(Abundance = Abundance.temp, 
                 SSB=SSB.temp, # need more work to correctly compute
                 Recruits=c(PopDy$Recruits,Stats.proj$R_proj),F_3d=F.temp, M_2d=M.temp) 
    
    # update model dimension and relevant quantities and generate new data sets (get ready for the subsequent assessment) 
    
    Year.proj.start = year_start + n_t
    Year.proj.end = Year.proj.start + n_proj -1 
    n_t = n_t + n_proj
    year_end = Year.proj.end
    
    # Call update files to simulate projected values
    # call input.data.update file
    source("C:...input.data.update.r")
    # call data.gen.proj file
    source("C:...data.gen.proj.r")
    
    # run AM
    unlink("NSLSAP01.STD") # remove the previous .std file in the folder
    runAD("NSLSAP01",verbose=T) # run AM
    # If the run converges, proceed
    fileNames=c("NSLSAP01.STD")
    if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
      file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
      do.AM = do.AM + 1
    }
    # If the run does not converge, record it, kick out that run and pick up wherever you left off
    while(file.exists(fileNames)==FALSE){ 
      noncon.counter = noncon.counter + 1
      ### Need to remove last three datapoints for nonconvergent runs
      survey.data <- survey.data[1:(n_t-n_proj),] # removes last 3 datapoints if no convergence (seems to work 8-17-2022)
      catch.data <- catch.data[-((3*n_t-n_proj+1):(3*n_t)),] # removes last three from fleet three
      catch.data <- catch.data[-((2*n_t-n_proj+1):(2*n_t)),] # should remove last three of fleet two
      catch.data <- catch.data[-((n_t-n_proj+1):n_t),]
      # call data.gen.proj again, and re-run the AM, repeat process until you get a convergent run
      source("C:...data.gen.proj.r")
      runAD("NSLSAP01",verbose=T)
      if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
        file.copy("NSLSAP01.REP", paste(scenario,'-',year_end,'-',iter,".","rep",sep=""),overwrite = TRUE) # copy converged .rep file and name it by scenario and iteration
        do.AM = do.AM + 1
      }
    }
    # Record updated "true" population dynamics
    obs.abun.at.size.array[1:(new_nt+n_proj),,do.AM+1,iter] <- PopDy$Abundance
    obs.ssb.array[1:(new_nt+n_proj),do.AM+1,iter] <- PopDy$SSB
    obs.rec.array[1:(new_nt+n_proj),do.AM+1,iter] <- PopDy$Recruits
  }
  # Check bias, as needed
  # report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
  # report <- read.rep('Fmsy-2020-1.rep') # or whatever the name of your scenario is
  # bias.check[iter] <- mean(apply(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100,2,mean))
  # bias.check[iter] <- max(abs(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100)) # maximum % bias of N-at-age
  iter = iter + 1 # update the iteration
} 

### LOOP TO RECORD MSE OUTPUTS ###

  iter = 1

  while(iter <= n.iter){ # must read in each report for each assessment over all iterations
    source("input.data.r") # get starting n_t, year_start values
    do.AM = 0
    while(do.AM < n.AM){
      report <- read.rep(paste(scenario,'-',year_end,'-',iter,".","rep",sep=""))
      # get F reference points from the AM 
      fyear <- report$BeginYear
      lyear <- report$EndYear
      new_nt <- length(fyear:lyear)
      
      F_fleets <- t(report$Fishing_Mortality[,new_nt])
      F_0.1 <- report$F0.1 # F0.1
      F_max <- report$Fmax # Fmax
      F_40spr <- report$F40SPR # F at 40% SPR, used as proxy for FMSY in this example
      F_msy <- F_spr40+f.err # add implementation error
      F_msy75 <- 0.75*F_msy 
      Msy_scen <- Msy_scen <- c(F_msy,F_msy75,F_0,F_0.1,F_max) 
      
      # update model dimensions
      Year.proj.start = year_start + n_t
      Year.proj.end = Year.proj.start + n_proj -1 
      n_t = n_t + n_proj
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
    
    # Check bias as needed
    # report <- read.rep('Fmsy-2020-1.rep') 
    bias.check[iter] <- mean(apply(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100,2,mean))
    # bias.check[iter] <- max(abs(((report$Abundance_at_Size-PopDy$Abundance[,,1])/PopDy$Abundance[,,1])*100)) # maximum % bias of N-at-age
    iter = iter + 1 # update the iteration
  }
  
  # calculate average annual variability for commercial fleet (not currently used, but it's there)
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
  
  abd.3d <- apply(abd.size.final, c(1,3), sum)
  exp.comm.final <- pred.catch.comm.final/abd.3d # commercial exploitation rates
  exp.rec.final <- pred.catch.rec.final/abd.3d # recreational exploitation rates, proxy for season length
  exp.disc.final <- pred.catch.disc.final/abd.3d
  catch.rate <- pred.catch.rec.final/exp.rec.final # catch rate
  # number of legal sized fish in the population
  abd.legal <- abd.size.final[,8:17,] # legal size classes for 11 inch size limit
  abd.legal.3d <- apply(abd.legal, c(1,3), sum)
  abd.legal.prop <- abd.legal.3d/abd.3d
  
  true.abd.size.final <- obs.abun.at.size.array[,,(n.AM+1),]
  true.abd.3d <- apply(true.abd.size.final, c(1,3), sum)
  
  # Convergence rate
  conv <- ((n.AM+1)*n.iter)/(((n.AM+1)*n.iter)+noncon.counter)
  
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
  
# These are just some quick and dirty diagnostic plots

# Jitter plots for population dynamics in the final assessment over all iterations
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
  }

  
# Jitter plots of catch by fleet
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

# Jitter plots of predicted vs. observed survey
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

# Abundance at size bubble plot
  abd.size.final <- abd.size.final[1:n_t,,] # need to remove final two years so uncertain recruitment doesn't mute size structure
  abd.size.final.mean <- rowMeans(abd.size.final, dim = 2) # takes the mean abundance at size over all iterations of the final assessment
  
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
  
  # Histograms for exploitation rate (season length), catch rate, and proportion legal
  windows(width=15,height=10)
  par(mfrow = c(1,3))
  hist(apply(exp.rec.final[1:n_t,1:n.iter],1,mean), xlab = "Mean exploitation rate (recreational)", main = "Exp rate/Season length")
  abline(v = median(apply(exp.rec.final[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)
  hist(apply(catch.rate[1:n_t,1:n.iter],1,mean), xlab = "Mean catch rate (recreational)", main = "Catch rate")
  abline(v = median(apply(catch.rate[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)
  hist(apply(abd.legal.prop[1:n_t,1:n.iter],1,mean), xlab = "Mean proportion of legal-sized fish (recreational)", main = "Proportion legal")
  abline(v = median(apply(abd.legal.prop[1:n_t, 1:n.iter],1,mean)), col = "red", lwd=3)