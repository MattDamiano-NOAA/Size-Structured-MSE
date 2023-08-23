# Size-structured MSE 
# Functions file
# M. Damiano
# Last updated: 8-23-2023

# Functions:

############## weight-length ###############
get_weight <- function(WL_pars_m,size_bm,n_t) 
{
  weight_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  for (t in 1:n_t)
  {
    weight_m[t,] = exp(WL_pars_m[t,1]+WL_pars_m[t,2]*log(size_bm)) 
  }
  return(weight_m)
}

############## maturity-sex change ###################
get_matu_prop <- function(fe_prop_pars_m,F50,size_bm,n_t) 
{
  maturity_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  per_fem = matrix(NA,nrow=n_t,ncol=length(size_bm))
  for (t in 1:n_t)
  {
    maturity_m[t,] = (1/(1+exp(-2*log(3)*(size_bm-fe_prop_pars_m[t,1])/fe_prop_pars_m[t,2])))
    # per_fem[t,] = 1-(1/(1+exp(-2*log(3)*(size_bm-F50)/fe_prop_pars_m[t,2]))) 
    # optional component if user wants to model exclusively female maturity, e.g., protogyny/dome-shaped maturity
  }
  maturity_m <- maturity_m # multiply maturity_m*per_fem
  return(maturity_m)
}

# optional fecundity function if assessment measures SSB in number of eggs

# get_fec <- function(fec_pars_m,size_bm,weight_m,n_t){
#   fec_mat = matrix(NA, nrow = n_t, ncol = length(size_bm))
#   for(t in 1:n_t){
#     fec_mat[t,] = exp(fec_pars_m[t,1]+fec_pars_m[t,2]*log(weight_m[t,]))
#   }
#   fec_mat <- fec_mat/weight_m # added to workaround SSB calcs issue / standardize fec
#   return(fec_mat)
# }

############## set natural mortality ############
get_M <- function(M_size_v,M_year_v) 
{
  
  n_l=length(M_size_v)
  n_t=length(M_year_v)
  
  M_m = array(NA, c(n_t,n_l,n_s))
  for (n in 1:n_t)
  {
    M_m[n,,] = M_year_v[n]*M_size_v 
  }
  M_m
  return(M_m) 
}

############# selectivities ###############
cal_sels <- function(sel_pars_v,sel_switch,dome_pars,size_bm) 
{
    if (sel_switch==2&length(sel_pars_v)!=2) 
    stop("check # of parameters")
  
    if (sel_switch==3&length(dome_pars)!=4)
    stop("check # of parameters for dome")
  
    sels <- array(NA, n_l)
    
    if (sel_switch==2) { # 2 if logistic/flat-topped, use logistic model
      sels <- 1.0/(1.0+exp((sel_pars_v[1]- size_bm)* sel_pars_v[2])); 
      sels <- sels/max(sels)
    }
    
    else { # 3 if dome-shaped, use double logistic model 
      sels_temp1 <- 1.0/(1.0+exp((dome_pars[1] - size_bm)* dome_pars[2]))
      sels_temp2 <- 1.0 - 1.0/(1.0+exp((dome_pars[3] - size_bm)* dome_pars[4]))
      sels <- (sels_temp1*sels_temp2)/max(sels_temp1*sels_temp2)
  } 
  return(sels)
}

###############
get_sels <- function(sel_pars_v,sel_switch,switch_v,n_fleets,size_bm,n_s,n_t,n_l)
{
  sels_4d = array(NA,c(n_t,n_l,n_s,n_fleets))
  for (f in 1:n_fleets){
    for (s in 1:n_s){
      for (t in 1:n_t){ 
        
          sels_4d[t,,s,f] = cal_sels(sel_pars_v,switch_v[f],dome_pars,size_bm) 
      }
    }
  }
  return(sels_4d) 
}

############## fishing mortality ##############
get_F <- function(F_array,sel_pars_v,sel_switch_v,n_fleets,size_bm,n_s,n_t,n_l)
{
  sels <- get_sels(sel_pars_v,sel_switch,switch_v,n_fleets,size_bm,n_s,n_t,n_l)
  F_4d = array(NA,c(n_t,n_l,n_s,n_fleets))
  
  for (f in 1:n_fleets){
    for (s in 1:n_s){
      for (t in 1:n_t){
        F_4d[t,,s,f] = F_array[t,s,f]*sels[t,,s,f] 
      }
    }
  }
  
  F_3d = array(0,c(n_t,n_l,n_s))
  temp = array(0,c(n_t,n_l,n_s))

  for (f in 1:n_fleets){ 
    if (n_s == 1){ 
      temp[,,1] = F_4d[,,,f]
    }else{
      temp = F_4d[,,,f]
    }
    F_3d = temp + F_3d
  }
  f_obj = list(F_3d=F_3d,F_4d=F_4d)
  return(f_obj)
}

get_F_inits <- function(F_array) 
{
  f_first_year_index      = matrix(0,nrow=n_fleets,ncol=n_s)       
  log_f_first_year        = matrix(0,nrow=n_fleets,ncol=n_s)
  log_F_devs              = array(0,dim=dim(F_array))
  
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    { 
      for (t in 1:n_t)
      {
        if (F_array[t,s,f]>0)
        {
          f_first_year_index[f,s] = t
          break
        }
      }
    }
  }
  f_first_year_index
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    {
      year.index.temp        = f_first_year_index[f,s]
      
      if (year.index.temp > 0){
        log_f_first_year[f,s]  = log(F_array[year.index.temp,s,f])
        
        for (t in 1:n_t)
        {
          if (t==f_first_year_index[f,s])
          {
            log_F_devs[t,s,f] = 0
          }
          if (t>f_first_year_index[f,s])
          {
            # log_F_devs[t,s,f] = F_devs_array[t,s,f] # If conditioned values are provided
            log_F_devs[t,s,f] = log(F_array[t,s,f])-log(F_array[t-1,s,f]) # if no conditioned values are provided
          }
        }
        break
      }
    }
  }
  f_inits.obj = list(log_f_first_year=log_f_first_year,log_F_devs=log_F_devs)
  return(f_inits.obj)
}

############## growth transition matrix #############
cal_GM <- function(Lmin,Lmax,Linc,Linf,SELinf,K,SEK,rhoLinfK)
{
  growmat <- NULL
  rhoLinfK <- abs(rhoLinfK)
  Ln <- seq(Lmin, Lmax+Linc, Linc)
  COV <- rhoLinfK * SELinf * SEK
  DL <- (Linf - Ln) * (1 - exp(-K))
  DL <- ifelse(DL < 0, 0, DL)
  VL <- SELinf^2 * (1 - exp(-K))^2 + (Linf - Ln)^2 * SEK^2 * 
    exp(-K)*exp(-K) - 2 * COV * (1 - exp(-K)) * (Linf - Ln) * 
    exp(-K)
  sqrt(VL)
  Ln+DL
  growmat <- matrix(0, nrow = length(Ln) - 1, ncol = length(Ln) - 1)
  
  for (L in 1:as.numeric(length(Ln) - 1)) {
    
    for (m in L:as.numeric(length(Ln) - 1)) {
      growmat[L, m] <- pnorm(Ln[m + 1]-(0.5*Linc), mean = Ln[L]+  
                               DL[L], sqrt(VL[m])) - pnorm(Ln[m]-(0.5*Linc), mean = Ln[L]+ 
                                                             DL[L], sqrt(VL[m]))
    }
  }
  growmat[is.na(growmat)] <- 0 # If there are any NAs, replace with 0s
  
  growmat <- growmat/rowSums(growmat, na.rm = T) 
  n.dim=(Lmax-Lmin)+Linc
  growmat_ful = matrix(0,nrow=n.dim,ncol=n.dim)
  
  if (nrow(growmat)<n.dim){
    growmat_ful[1:nrow(growmat),1:ncol(growmat)]=growmat
    for(nn in 1:(n.dim-nrow(growmat))){
      growmat_ful[nrow(growmat)+nn,ncol(growmat)+nn]=1
    }
  }
  else{
    growmat_ful=growmat[1:n.dim,1:n.dim]
  }
  return(growmat)
} 

get_GM <- function(n_l,n_growblock,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
{
  grow_3d = array(0,c(n_l,n_l,n_growblock))
  grow_3d
  for (t in 1:n_growblock){
    grow_3d[,,t] = cal_GM(Lmin,Lmax,Linc,Linf_v[t],SELinf_v[t],K_v[t],SEK_v[t],rhoLinfK_v[t])
  }
  return(grow_3d)
}

####### Stock-recruitment relationship ######
get_R <- function(SSB,alpha,beta,indicator)
{
  if (indicator==1)R=alpha # no functional relationship, aka "mean recruitment model"
  if (indicator==2)R=alpha*SSB/(beta+SSB) # Beverton-Holt model
  if (indicator==3)R=alpha*SSB*exp(-beta*SSB) # Ricker model
  return(R)
}

get_ssb.fracs <- function(ssb_month,n_s) #This creates periods during the year when spawning biomass is contributed
{
  SSB_frac_v    = c()
  SSB_frac_v[1] = (ssb_month-1)/12
  temp=SSB_frac_v[1] *12;temp1=0
  for(i in 1:n_s)
  { 
    temp1=temp1+12/n_s
    if(temp<=temp1)
    {
      SSB_frac_v[2] = i
      SSB_frac_v[3] = (temp-(temp1-12/n_s))/(12/n_s)
      break
    }
  }
  return(SSB_frac_v)
}

############### pop dyn ##########################
pop_dyn <- function(N_init_tot,N_init_comp,F_3d,M_m,SSB_frac_v,maturity_m,weight_m,alpha,beta,indicator,R_devs_v,R_l_pro,n_t,n_s,n_growblock,growth_array)
{ 
  N_3d = array(0,c(n_t,n_l,n_s))
  SSB  = c()
  R    = c()
  
  for (t in 1:n_t){
    
    if(t>1)
      R[t] = get_R(SSB[t-1],alpha,beta,indicator)*exp(R_devs_v[t])
    
    for (s in 1:n_s){
      
      if(t==1&s==1){
        
        # Set up values for the first year (and season if applicable)
        
        # If there is seasonality in the model, use these to initialize F and M
        # F_temp_year = apply(F_3d[1,,1:n_s],1,sum)
        # M_temp_year = apply(M_m[1,,1:n_s],1,sum)
        
        # Mortality
        F_temp_year = F_3d[1,,1] 
        M_temp_year = M_m[1,,] 
        
        # Initial abundance at size
        N_3d[1,,1] = N_init_tot*N_init_comp #Starting N numbers by length
        # Initial SSB
        SSB[1] = sum(N_3d[1,,1]*exp(-SSB_frac_v[3]*(F_temp_year+M_temp_year))*maturity_m[1,]*weight_m[1,]) # fecundity can be added here to calculate SSB in eggs
        # Initial recruitment
        R[1] = get_R(SSB[1],alpha,beta,indicator)*exp(R_devs_v[1])
        # Add recruitment at size to initial abundance at size
        N_3d[1,,1] = N_3d[1,,1] + R[1]*R_l_pro
      }
      # for each year in the first season (also the main N equation if there's no seasonality)
      if(t>1&s==1){
        N_3d[t,,s]=(N_3d[t-1,,n_s]*exp(-(F_3d[t-1,,n_s]+M_m[t-1,,1])))%*%growth_array[,,t-1]+R[t]*R_l_pro #%*% operator is for matrix mult
        n_growblock[t-1] # not currently in use
        if(SSB_frac_v[2]==1&n_s>1){ # only used if there is seasonality
          SSB[t] = sum(N_3d[t,,1]*exp(-SSB_frac_v[3]*(F_3d[t,,]+M_m[t,]))*maturity_m[t,]*weight_m[t,])
        }
      }
      
      if(n_s==1&t>1){ # default SSB equation
        SSB[t] = sum(N_3d[t,,n_s]*exp(-SSB_frac_v[1]*(F_3d[t,,n_s]+M_m[t,,1]))*maturity_m[t,]*weight_m[t,])
      }
    }
  }
  #If there is seasonality in the model configuration 
  if(s>1){
    N_3d[t,,s]=(N_3d[t,,s-1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1])))%*%grow_3d[,,growblock]+R[t]*R_l_pro
    
    if(s==SSB_frac_v[2]&n_s>1){
      SSB[t] = sum(N_3d[t,,s]*exp(-SSB_frac_v[3]*(F_3d[t,,s]+M_m[t,,s]))*maturity_m[t,]*weight_m[t,])
    }
  }
  pop.dyn.obj = list(Abundance=N_3d,SSB=SSB,Recruits=R,F_3d=F_3d,M_2d=M_m)
  return(pop.dyn.obj)
}


############### calulate catch ##########################
cal_catch <- function(f_obj,M_m,N_3d,n_t,n_l,n_s,n_fleets)
{
  catch_4d  = array(NA,c(n_t,n_l,n_s,n_fleets))
  catch_3d  = array(NA,c(n_t,n_fleets,n_s))
  catch_tot = matrix(NA,nrow=n_t,ncol=n_fleets)
  
  F_3d = f_obj$F_3d
  F_4d = f_obj$F_4d
  
  for (f in 1:n_fleets)
  {
    for (t in 1:n_t)
    {
      for (s in 1:n_s)
      {
        catch_4d[t,,s,f] = N_3d[t,,s]*(1-exp(-(F_3d[t,,s]+M_m[t,,s])))*(F_4d[t,,s,f]/(F_3d[t,,s]+M_m[t,,s]))#Just the baranov catch equation, easy
        catch_3d[t,f,s] = sum(catch_4d[t,,s,f])
      }
    }
  }
  catch.obj = list(catch_4d=catch_4d,catch_3d=catch_3d)
  return(catch.obj)   
}

############### generate survey data (no error) ##########################
get_survey.fracs <- function(survey_months,n_survey,n_s) # calculates when in the year the survey occurs
{
  survey_frac_m    = matrix(NA,nrow=n_survey,ncol=3) 
  for (n in 1:n_survey)
  {
    survey_frac_m[n,1] = survey_months[n]/12
    
    temp=(survey_months[n]/12)*12
    temp1=0
    
    for(i in 1:n_s)
    {
      temp1=temp1+12/n_s
      if (temp<=temp1)
      {
        survey_frac_m[n,2] = i
        survey_frac_m[n,3] = (temp-(temp1-12/n_s))/(12/n_s)
        break
      }
    }
    
  }
  return(survey_frac_m)
}

# calculates simulated indices of relative abundance
cal_survey_inds <- function(N_3d,F_3d,M_m,survey_frac_m,survey_q,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,n_t,n_l)
{  
  survey.temp  = array(0,c(n_t,n_l,n_survey))
  survey.sels  = matrix(NA,nrow=n_survey,ncol=n_l)
  survey.obj   = list()
  year_start = survey_periods[,1]
  year_end = survey_periods[,2]
  year_start
  year_end
  
  for (n in 1:n_survey) # if there is more than one survey
  {
    survey.sels[n,] = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)
    n1 = survey_periods[n,1] - year_start + 1
    n2 = n_t - (year_end-year_start)+1

    for (t in 1:n_t)
    {
      if(n_s <= 1) 
      {
        survey.temp[t,,n] = N_3d[t,,1]*exp(-survey_frac_m[n,1]*(F_3d[t,,1]+M_m[t,,1]))*survey_q[n]*survey.sels[n,]
      }
      if(n_s > 1)
      {
        temp = survey.fracs[n,2]
        survey.temp[t,,n] = N_3d[t,,temp]*exp(-survey_frac_m[n,3]*(F_3d[t,,temp]+M_m[t,,temp]))*survey_q[n]*survey.sels[n,]
      }
      
    }
    survey.obj[[n]] = cbind(apply(survey.temp[,,n],1,sum),survey.temp[,,n]) 
  }
  return(survey.obj)
  
}

############### generate data (with random error) ##########################
add.error.lognormal <- function(data,cv)
{
  sd               = sqrt(log((cv^2+1)))
  data.error       = rlnorm(length(data),log(data),sd)
  return(data.error)
}

add.error.multino <- function(data,ess)
{
  dim.data      = dim(data)  
  data.error    = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
  sample.temp   = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
  
  original.com=sweep(data,1,rowSums(data),"/")
  
  for (i in 1:dim.data[1])
  {
    if (sum(data[i,])>0)
    {
      sample.temp[i,] = rmultinom(1, ess[i], prob=original.com[i,])
      data.error[i,]  = sample.temp[i,]/sum(sample.temp[i,])
    }  
  }
  return(data.error)
}

# Generate catch data in format required for the size-structured assessment model
catch.data.gen <- function(year_start,year_end,catch.obj,catch.cv,catch.ess,add.error)
{
  
  n.season       = lapply(catch.obj, dim)$catch_4d[3]
  n.fleets       = lapply(catch.obj, dim)$catch_4d[4]
  n.lengthbin    = lapply(catch.obj, dim)$catch_4d[2]
  n.year         = lapply(catch.obj, dim)$catch_4d[1]
  
  n_s = n.season
  n_t = n.year
  n_fleets = n.fleets
  
  year_start = year_start
  year_end = year_end
  catch.data.m   = matrix(-1,nrow=n_t*n_s*n_fleets,ncol=9+length(size_bm)) #changed n.lengthbins to what I've been using to count them 'size_bm'
  catch.data.m
  length(size_bm)
  
  year      = rep(seq(year_start,year_end,1),n_s*n_fleets)
  season    = rep(rep(1:n_s,each=n_t),n_fleets)
  fleet     = rep(rep(1:n_fleets,each=n_t*n_s))
  catch.tot = c()
  cv        = c()
  cpue      = rep(-1,n_t*n_s*n_fleets)
  cpue.flag = rep(1,n_t*n_s*n_fleets)
  cpue.cv   = rep(-1,n_t*n_s*n_fleets)
  ess       = c()
  comp      = data.frame()
  
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    { 
      cv.temp          = catch.cv[,f,s]
      cv               = c(cv,cv.temp)
      catch.tot.temp   = catch.obj$catch_3d[,f,s]   
      catch.tot        = c(catch.tot,catch.tot.temp)
      ess.temp         = catch.ess[,f,s]
      ess              = c(ess,ess.temp)
      comp.temp        = catch.obj$catch_4d[,,s,f]
      comp             = rbind(comp,comp.temp)
      
      if(add.error==TRUE)
      {
        catch.tot    = add.error.lognormal(catch.tot,cv)
        comp         = add.error.multino(comp,ess)
      }
      
    }
  }
  
  catch.data.m[,1]                   = year
  catch.data.m[,2]                   = season
  catch.data.m[,3]                   = fleet
  catch.data.m[,4]                   = catch.tot
  catch.data.m[,5]                   = cv
  catch.data.m[,6]                   = cpue
  catch.data.m[,7]                   = cpue.flag
  catch.data.m[,8]                   = cpue.cv
  catch.data.m[,9]                   = ess
  catch.data.m[,10:(9+n.lengthbin)]  = as.matrix(comp)
  return(catch.data.m)
}

# Generate survey index data in the format required for the size-structured assessment model
survey.data.gen <- function(survey_periods,survey_months,n_survey,survey.obj,survey.cv,survey.ess,add.error)
{
  nrow           = 0
  nrow.temp      = 0
  year           = c()
  indexN         = c()
  month          = c()
  index          = c()
  cv             = c()
  ess            = c()
  comp           = data.frame()
  n_sizeb        = lapply(survey.obj, dim)[[1]][2]
  
  for (i in 1:n_survey)
  {
    nrow.temp    = lapply(survey.obj, dim)[[i]][1]
    nrow         = nrow + nrow.temp
    year.temp    = seq(survey_periods[i,1],survey_periods[i,2],1)
    year         = c(year,year.temp)
    indexN.temp  = rep(i,nrow.temp)
    indexN       = c(indexN,indexN.temp)
    month.temp   = rep(survey_months[i],nrow.temp)
    month        = c(month,month.temp)
    index.temp   = survey.obj[[i]][,1]
    index        = c(index,index.temp)
    cv.temp      = survey.cv[[i]]
    cv           = c(cv,cv.temp)
    ess.temp     = survey.ess[,i,]
    ess          = c(ess,ess.temp)
    comp.temp    = survey.obj[[i]][,2:n_sizeb]
    comp         = rbind(comp,comp.temp)
    
    if(add.error==TRUE)
    {
      index      = add.error.lognormal(index,cv)
      comp       = add.error.multino(comp,ess)
    } 
    
  }
  survey.data.m  = matrix(-1,nrow=length(index),ncol=(6+n_sizeb-1))
  
  survey.data.m[,1]                 = year
  survey.data.m[,2]                 = indexN
  survey.data.m[,3]                 = month
  survey.data.m[,4]                 = index
  survey.data.m[,5]                 = cv
  survey.data.m[,6]                 = ess
  survey.data.m[,7:(6+n_sizeb-1)]   = as.matrix(comp)
  
  return(survey.data.m)  
}