# Dolphinfish MSE
# Dolphinfish OM: functions file
# Last modified: 10-12-2023 by M. Damiano

# This file calls all of the functions outlined in the technical documentation
# See the input.data file for descriptions of parameters fed to each function

############## weight-length ###############
# Calculates the weight-length relationship that is used later for spawning stock biomass (SSB)
get_weight <- function(WL_pars_m,size_bm,n_t)
{
  weight_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  for (t in 1:n_t)
  {
    weight_m[t,] = exp(WL_pars_m[t,1] + WL_pars_m[t,2]*log(size_bm)) #should create a year x size matrix of weight at length (k) for each year, time invariant
  }
  return(weight_m)
}

############## maturity ###################
# Calculates the number of mature fish - this is typically females, which
# is why it uses the matrix of proportion female (fe_prop_pars_m)
get_matu_prop <- function(fe_prop_pars_m,size_bm,n_t) #fe_prop_pars_m is a 1x2 matrix of the 2 maturity parameters
{
  maturity_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  # per_fem = matrix(NA,nrow=n_t,ncol=length(size_bm)) # remnant from using this for sex change in protogynous fish
  for (t in 1:n_t)
  {
    maturity_m[t,] = (1 / (1 + exp(-2 * log(3) * (size_bm - fe_prop_pars_m[t,1]) / fe_prop_pars_m[t,2])))
    # per_fem[t,] = 1-(1/(1+exp(-2*log(3)*(size_bm-F50)/fe_prop_pars_m[t,2])))
  }
  return(maturity_m)
}

# Optional fecundity function if fecundity information is available

# get_fec <- function(fec_pars_m,size_bm,weight_m,n_t){
#   fec_mat = matrix(NA, nrow = n_t, ncol = length(size_bm))
#   for(t in 1:n_t){
#     fec_mat[t,] = exp(fec_pars_m[t,1]+fec_pars_m[t,2]*log(weight_m[t,]))
#   }
#   fec_mat <- fec_mat/weight_m # added to workaround SSB calcs issue / standardize fec
#   return(fec_mat)
# }

############## set natural mortality ############
# M function can technically scale with age or size
get_M <- function(M_size_v,M_year_v) #vectors of mortality by season, size and year
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
# The selectivity can either be logistic or dome-shaped (double logistic)
# Uses the switch to determine if the right number of parameters have been input
cal_sels <- function(sel_pars_v,sel_switch,dome_pars,size_bm) #pars_v contains the two parameters for a logistic selectivity curve. This is fine given that fleets for BSB follow this pattern
{
  if (sel_switch==2&length(sel_pars_v)!=2)
    stop("check # of parameters")

  if (sel_switch==3&length(dome_pars)!=4)
    stop("check # of parameters for dome")

  sels <- array(NA, n_l)

  if (sel_switch==2) { # if logistic/flat-topped
    sels <- 1.0/(1.0+exp((sel_pars_v[1]- size_bm)* sel_pars_v[2])); #Gives a selectivity by length, changes pars_v[2] to negative to reflect logistic curve.
    #Think of slope as how fast sel increases with size
    sels <- sels/max(sels)
  }

  else { # if dome-shaped, use double logistic model (Method 1990)
    sels_temp1 <- 1.0/(1.0+exp((dome_pars[1] - size_bm)* dome_pars[2]))
    sels_temp2 <- 1.0 - 1.0/(1.0+exp((dome_pars[3] - size_bm)* dome_pars[4]))
    sels <- (sels_temp1*sels_temp2)/max(sels_temp1*sels_temp2)
  }
  sels
  return(sels)
}

###############

# This function used to calculate an array of selectivities to use
# in the calculation of fishing mortality, but because we have so many
# different selectivity patterns, the general function above cannot meet
# those needs. So, the sels_4d array is assembled in the input.data file
# and treated as a parameter for the fishing mortality function.

# get_sels <- function(sel_pars_v,sel_switch,switch_v,n_fleets,size_bm,n_s,n_t,n_l)
#   #create a selectivity array of n_t x n_l matrices by fleet
# {
#   # dim pars_m - c(n_block,6)
#   #block_array <- array(NA, c(n_t,n_s,n_block))
#   #block_array
#   # swith_v swith for each block
#   sels_4d = array(NA,c(n_t,n_l,n_s,n_fleets))
#   #sels_4d
#   for (f in 1:n_fleets){ # consider areas as fleets approach to different selectivities in the general private rec fleet
#     for (s in 1:n_s){
#       for (t in 1:n_t){ # needs the block back so it can switch between sel
#
#           sels_4d[t,,s,f] = cal_sels(sel_pars_v,switch_v[f],dome_pars,size_bm) #now it just calls the same selectivity parameters (logistic, 2 param)
#       }
#     }
#   }
#   sels_4d
#   return(sels_4d)
# }

############## fishing mortality ##############
# Calculates fishing mortality at size (F*selectivity) for each fleet
# Can be annual or seasonal
get_F <- function(F_array,sels_4d,n_fleets,size_bm,n_s,n_t,n_l)
{
  # dim(F_m) - c(n_t,n_s,n_fleets)
  sels <- sels_4d

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
  # for (f in 1:n_fleets){
  #   temp = F_4d[,,,f]
  #   F_3d = temp + F_3d
  # }

  for (f in 1:n_fleets){
    if (n_s == 1){ # Whether 4D or 3D depending on if season is set to 1 for annual time steps
      temp[,,1] = F_4d[,,,f]
    }else{
      temp = F_4d[,,,f]
    }
    F_3d = temp + F_3d
  }

  f_obj = list(F_3d=F_3d,F_4d=F_4d)
  return(f_obj)
}

# This primary use of this function is to get year to year or season to season
# deviations in F; the rest is a remnant of the old model that needed it
# for the estimation components

get_F_inits <- function(F_array)
{
  f_first_year_index      = matrix(0,nrow=n_fleets,ncol=n_s)       # dim - c(n_fleets,season)
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
            # log_F_devs[t,s,f] = F_devs_array[t,s,f] # For conditioned values provided
            log_F_devs[t,s,f] = log(F_array[t,s,f])-log(F_array[t-1,s,f]) # if no conditioned values provided
          }
        }
        break
      }
    }
  }
  log_F_devs
  f_inits.obj = list(log_f_first_year=log_f_first_year,log_F_devs=log_F_devs)
  f_inits.obj
  return(f_inits.obj)

}

############## growth transition matrix #############
# Calculates the expected growth increment for each size bin based on
# Von bertalanffy growth function parameters, and generates a k x k matrix
# of growth probabilities - cornerstone of size-based assessments
cal_GM <- function(Lmin,Lmax,Linc,Linf,SELinf,K,SEK,rhoLinfK)
{
  growmat <- NULL
  rhoLinfK <- abs(rhoLinfK)
  Ln <- seq(Lmin, Lmax+Linc, Linc)
  COV <- rhoLinfK * SELinf * SEK
  #Delta L, change in length
  #The growth increment
  DL <- (Linf - Ln) * (1 - exp(-K))
  DL <- ifelse(DL < 0, 0, DL)
  #Variance of delta L, 11-1-19 changed - 2*COV to + 2*COV per equation in 2017 paper
  # 6/26/2020: changed it back to be consistent with EM
  VL <- SELinf^2 * (1 - exp(-K))^2 + (Linf - Ln)^2 * SEK^2 *
    exp(-K)*exp(-K) - 2 * COV * (1 - exp(-K)) * (Linf - Ln) *
    exp(-K)
  VL
  sqrt(VL)
  Ln+DL
  #7-9-2020: math is consistent with EM for calculating variance
  growmat <- matrix(0, nrow = length(Ln) - 1, ncol = length(Ln) - 1)
  growmat
  for (L in 1:as.numeric(length(Ln) - 1)) {

    for (m in L:as.numeric(length(Ln) - 1)) {
      growmat[L, m] <- pnorm(Ln[m + 1]-(0.5*Linc), mean = Ln[L]+  # add -0.5*Linc to get midpoint size
                               DL[L], sqrt(VL[m])) - pnorm(Ln[m]-(0.5*Linc), mean = Ln[L]+ # add -0.5*Linc
                                                             DL[L], sqrt(VL[m]))
    }
  }
  growmat
  #11-1-19 need to replace row and column 20 with zeros b/c 0 probability of growing if they are as large as they can be
  growmat[is.na(growmat)] <- 0
  growmat
  #11-1-19 every row should sum to 1 by rules of probability

  growmat <- growmat/rowSums(growmat, na.rm = T)
  n.dim=(Lmax-Lmin)+Linc
  growmat_ful = matrix(0,nrow=n.dim,ncol=n.dim)

  if (nrow(growmat)<n.dim){
    growmat_ful[1:nrow(growmat),1:ncol(growmat)]=growmat
    for(nn in 1:(n.dim-nrow(growmat))){
      growmat_ful[nrow(growmat)+nn,ncol(growmat)+nn]=1
    }
  }else{
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

############## Stock Recruitment ################
# Selects and calculates the stock recruitment relationship, if there is one
get_R <- function(SSB,alpha,beta,indicator)
{
  if (indicator==1)R=alpha #no functional relationship
  if (indicator==2)R=alpha*SSB/(beta+SSB) #BH
  if (indicator==3)R=alpha*SSB*exp(-beta*SSB) #Ricker
  return(R)
}
# This creates periods during the year when spawning biomass is contributed
get_ssb.fracs <- function(ssb_month,n_s)
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
# Central population dynamics function
# Calculates abundance (N) at size, recruitment (R), and SSB
collect_population_dynamics <- function(N_init_tot,N_init_comp,F_3d,M_m,SSB_frac_v,maturity_m,weight_m,alpha,beta,indicator,R_devs_v,R_l_pro,years_count,n_s,n_r,n_growblock,growth_array,mov.mat,mov.mat2,mov.mat3,mov.mat4)
{
  N_3d = array(0,c(years_count,n_l,n_s,n_r))
  SSB  = c()
  R    = c()
  # s=2 # for testing
  for (t in 1:years_count){

    if(t>1)
      R[t] = get_R(SSB[t-1],alpha,beta,indicator)*exp(R_devs_v[t])

    for (s in 1:n_s){

      if(t==1&s==1){
        # Values for the first time step
        F_temp_year = F_3d[1,,s]
        M_temp_year = M_m[1,,s]
        #Checks
        # M_m
        # F_temp_year
        # M_temp_year
        # F_4d
        exp(-SSB_frac_v[3]*(F_temp_year+M_temp_year))

        N_3d[1,,s,] = N_init_tot*N_init_comp # Starting N numbers by length, treat as same for all regions
        # N_3d
        #SSB[1] = sum(N_3d[1,,1]*exp(-SSB_frac_v[1]*(F_temp_year+M_temp_year))*PropFemale_m[1,]*weight_m[1,])
        SSB[1] = sum(N_3d[1,,s,]*exp(-SSB_frac_v[3]*(F_temp_year+M_temp_year))*maturity_m[1,]*weight_m[1,])
        # -- to match up with assessment model (SSB1 will only be correct when spawning time is within season1)
        SSB
        R[1] = get_R(SSB[1],alpha,beta,indicator)*exp(R_devs_v[1])
        R
        N_3d[1,,s,] = N_3d[1,,s,] + R[1]*R_l_pro
        N_3d[1,,1,]

      }
      # Values for all years in the first season
      if(t>1&s==1){
        # calculates abundance for each region by adding those that stay in the region to those that move to that region
        N_3d[t,,s,1]=(N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[1,2:7]))+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,2]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,3]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,4]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,6]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[1,7]

        N_3d[t,,s,2]=(N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[2,1],mov.mat[2,3:7]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,1]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,3]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,4]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,6]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[2,7]

        N_3d[t,,s,3]=(N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[3,1:2],mov.mat[3,4:7]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,1]+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,2]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,4]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,6]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[3,7]

        N_3d[t,,s,4]=(N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[4,1:3],mov.mat[4,5:7]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,1]+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,2]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,3]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,6]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[4,7]

        N_3d[t,,s,5]=(N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[5,1:4],mov.mat[5,6:7]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,1]+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,2]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,3]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,6]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[5,7]

        N_3d[t,,s,6]=(N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[6,1:5],mov.mat[6,7]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,1]+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,2]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,3]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,4]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,5]+
          (N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[6,7]

        N_3d[t,,s,7]=(N_3d[t-1,,s,7]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*(1-sum(mov.mat[7,1:6]))+
          (N_3d[t-1,,s,1]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,1]+
          (N_3d[t-1,,s,2]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,2]+
          (N_3d[t-1,,s,3]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,3]+
          (N_3d[t-1,,s,4]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,4]+
          (N_3d[t-1,,s,5]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,5]+
          (N_3d[t-1,,s,6]*exp(-(F_3d[t-1,,s]+M_m[t-1,,s]))%*%growth_array[,,t-1]+R[t]*R_l_pro)*mov.mat[7,6]

        # n_growblock[t-1]
        if(SSB_frac_v[2]==1&n_s>1){
          SSB[t] = sum(N_3d[t,,s,]*exp(-SSB_frac_v[3]*(F_3d[t,,s]+M_m[t,,s]))*maturity_m[t,]*weight_m[t,])
        }
      }

      if(n_s==1&t>1){
        SSB[t] = sum(N_3d[t,,s,]*exp(-SSB_frac_v[1]*(F_3d[t,,s]+M_m[t,,s]))*maturity_m[t,]*weight_m[t,])
      }

      #If there is seasonality in the model config.
      if(s==2){
        # s=2
        # t=1
        N_3d[t,,s,1]=(N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[1,2:7]))+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[1,7]

        N_3d[t,,s,2]=(N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[2,1],mov.mat2[2,3:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,1]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[2,7]

        N_3d[t,,s,3]=(N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[3,1:2],mov.mat2[3,4:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,2]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[3,7]

        N_3d[t,,s,4]=(N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[4,1:3],mov.mat2[4,5:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,3]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[4,7]

        N_3d[t,,s,5]=(N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[5,1:4],mov.mat2[5,6:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,4]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[5,7]

        N_3d[t,,s,6]=(N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[6,1:5],mov.mat2[6,7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,5]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[6,7]

        N_3d[t,,s,7]=(N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat2[7,1:6]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat2[7,6]
      }
      if(s==3){
        # s=2
        # t=1
        N_3d[t,,s,1]=(N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[1,2:7]))+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[1,7]

        N_3d[t,,s,2]=(N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[2,1],mov.mat2[2,3:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,1]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[2,7]

        N_3d[t,,s,3]=(N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[3,1:2],mov.mat2[3,4:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,2]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[3,7]

        N_3d[t,,s,4]=(N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[4,1:3],mov.mat2[4,5:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,3]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[4,7]

        N_3d[t,,s,5]=(N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[5,1:4],mov.mat2[5,6:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,4]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[5,7]

        N_3d[t,,s,6]=(N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[6,1:5],mov.mat2[6,7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,5]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[6,7]

        N_3d[t,,s,7]=(N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat3[7,1:6]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat3[7,6]
      }
      if(s==4){
        # s=2
        # t=1
        N_3d[t,,s,1]=(N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[1,2:7]))+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[1,7]

        N_3d[t,,s,2]=(N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[2,1],mov.mat4[2,3:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,1]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[2,7]

        N_3d[t,,s,3]=(N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[3,1:2],mov.mat4[3,4:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,2]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[3,7]

        N_3d[t,,s,4]=(N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[4,1:3],mov.mat4[4,5:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,3]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[4,7]

        N_3d[t,,s,5]=(N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[5,1:4],mov.mat4[5,6:7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,4]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,6]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[5,7]

        N_3d[t,,s,6]=(N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[6,1:5],mov.mat4[6,7]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,5]+
          (N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[6,7]

        N_3d[t,,s,7]=(N_3d[t,,s-1,7]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*(1-sum(mov.mat4[7,1:6]))+
          (N_3d[t,,s-1,1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,1]+
          (N_3d[t,,s-1,2]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,2]+
          (N_3d[t,,s-1,3]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,3]+
          (N_3d[t,,s-1,4]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,4]+
          (N_3d[t,,s-1,5]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,5]+
          (N_3d[t,,s-1,6]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1]))%*%growth_array[,,t]+R[t]*R_l_pro)*mov.mat4[7,6]
      }
      if(s==SSB_frac_v[2]&n_s>1){ # for now, SSB is summed across all regions and season to an annual metric
        SSB[t] = sum(N_3d[t,,s,]*exp(-SSB_frac_v[3]*(F_3d[t,,s]+M_m[t,,s]))*maturity_m[t,]*weight_m[t,])
      }
    }
  }
  pop.dyn.obj = list(Abundance=N_3d,SSB=SSB,Recruits=R,F_3d=F_3d,M_2d=M_m)
  return(pop.dyn.obj)
}

############### calulate catch ##########################
# Calculates catch at size for each region
cal_catch <- function(f_obj,M_m,N_3d,years_count,n_l,n_s,n_fleets,n_r)
{
  catch_4d  = array(NA,c(years_count,n_l,n_s,n_fleets,n_r))
  catch_3d  = array(NA,c(years_count,n_s,n_fleets))
  catch_tot = matrix(NA,nrow=years_count,ncol=n_fleets)

  F_3d = f_obj$F_3d
  F_4d = f_obj$F_4d

  for(r in 1:n_r){
    for (f in 1:n_fleets)
    {
      for (t in 1:years_count)
      {
        for (s in 1:n_s)
        {
          catch_4d[t,,s,f,r] = N_3d[t,,s,r]*(1-exp(-(F_3d[t,,s]+M_m[t,,s])))*(F_4d[t,,s,f]/(F_3d[t,,s]+M_m[t,,s]))
          catch_3d[t,s,f] = sum(catch_4d[t,,s,f,])
        }
      }
    }
  }
  catch.obj = list(catch_4d=catch_4d,catch_3d=catch_3d)
  return(catch.obj)
}

############### generate fishery-dependent CPUE data ##########################
# This uses the generic effort time series from the input and catch time series
# from the last function to calculate catch per unit effort (CPUE)
get_cpue <- function(catch_3d, eff_switch_v, eff_array){
  cpue_array <- array(NA, dim = c(years_count,n_s,n_fleets))
  for(i in 1:length(eff_switch_v)){
    if(eff_switch_v[i]==1){
      cpue_array[,,i] <- catch_3d[,,i]/eff_array[,1]
    }
    if(eff_switch_v[i]==2){
      cpue_array[,,i] <- catch_3d[,,i]/eff_array[,2]
    }
    if(eff_switch_v[i]==3){
      cpue_array[,,i] <- catch_3d[,,i]/eff_array[,3]
    }
    if(eff_switch_v[i]==4){
      cpue_array[,,i] <- catch_3d[,,i]/eff_array[,4]
    }
    if(eff_switch_v[i]==5){
      cpue_array[,,i] <- catch_3d[,,i]/eff_array[,5]
    }
  }
  return(cpue_array)
}

# This function has been commented out because we needed some more bespoke
# functionality to fit OM-generated data to observed indices of relative abundance

############### generate survey data (no error) ##########################
# this will need to change to reflect when "indices" from VAST are generated, i.e., by season. Maybe it's enough to change this to
# get_survey.fracs <- function(survey_months_v,n_survey,n_s)
# {
#   survey_frac_v    = c()
#   # for (n in 1:n_survey)
#   # {
#     # This newer code sets up quarterly fractions for applying survival in the cal_survey_inds function below
#     survey_frac_v[1] = (survey_months_v[1]+1)/12
#     survey_frac_v[2] = (survey_months_v[2]+1)/12
#     survey_frac_v[3] = (survey_months_v[3]+1)/12
#     survey_frac_v[4] = (survey_months_v[4]+1)/12
#     # temp=(survey_months[n]/12)*12
#     # temp1=0
#
#     # for(i in 1:n_s)
#     # {
#     #   temp1=temp1+12/n_s
#     #   if (temp<=temp1)
#     #   {
#     #     survey_frac_m[n,2] = i
#     #     survey_frac_m[n,3] = (temp-(temp1-12/n_s))/(12/n_s)
#     #     break
#     #   }
#     # }
#
#   # }
#   return(survey_frac_v)
# }

# Calculates "survey" index, i.e., index of relative abundance, based on
# N, F, M, and assumed "survey" selectivity for each region
cal_survey_inds <- function(N_3d,F_3d,M_m,survey_frac_v,survey_q_v,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,years_count,n_l,n_r)
{
  survey.temp  = array(0,c(years_count,n_l,n_survey,n_r)) # needs spatial dimension and season dimension. I think this whole function needs to be re-written
  survey.sels  = matrix(NA,nrow=n_survey,ncol=n_l)
  survey.obj   = array(0,c(years_count,n_l+1,n_survey,n_r))
  # year_start = survey_periods[,1]
  # year_end = survey_periods[,2]
  # year_start
  # year_end

  for (n in 1:n_survey){
    for (r in 1:n_r){
      survey.sels[n,] = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm) # could probably assume the same selectivity but different q for VAST indices
      # survey.sels #looks good
      # n1 = survey_periods[n,1] - year_start + 1
      # n2 = years_count - (year_end-year_start)+1
      # n1
      # n2

      for (t in 1:years_count)
      {
        for(s in 1:n_s){
          survey.temp[t,,n,r] = N_3d[t,,n,r]*exp(-survey_frac_v[s]*(F_3d[t,,s]+M_m[t,,s]))*survey_q_v[r]*survey.sels[n,]
          # fitting goes here?
        }
      }
      survey.obj[,,n,r] = apply(survey.temp[,,n,r],1,sum)
    }
    #need to come up with a different storage object
  } # also we are not fitting to index length comps, so we can get rid of that
  return(survey.obj)
}

# fitting function

fit.index <- function(pars,ind.obs,N_3d,F_3d,M_m,survey_frac_v,survey_sel_switch,survey_sel_pars,dome_pars,size_bm,years_count,n_l,s,r)
{
  pars[1] -> q
  ind.temp = array(NA,dim=c(years_count,n_l))
  ind.pred <- c()
  survey_sels = matrix(NA,ncol=1,nrow=n_l)
  survey.sels = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)
  for (t in 1:years_count){
    ind.temp[t,] <- N_3d[t,,s,r]*exp(-survey_frac_v[s]*(F_3d[t,,s]+M_m[t,,s]))*q*survey.sels # if we're using the N_3d object, then they already survived...
  }
  ind.pred <- apply(ind.temp,1,sum)
  obj = sum((ind.obs-ind.pred)^2)
  return(obj)
}

# Function that does the actual fitting to observed indices/CPUE
# Currently uses sum of squares to estimate a single catchability "q"
# for each region and each season (total = 28)
get.fitted.pred.index <- function(ind.vast,N_3d,F_3d,M_m,survey_frac_v,survey_q_v,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,years_count,n_l,n_r){
  q.vec <- array(NA,dim=c(n_s,n_r))
  ind.pred.temp <- array(NA,dim=c(years_count,n_l,n_s,n_r)) # use size-structured dynamics and selectivity to generate predicted index
  ind.pred <- array(NA,dim=c(years_count,n_s,n_r)) # there's no length info in the index, so output needs to be aggregated
  survey_sels = matrix(NA,ncol=1,nrow=n_l)
  survey.sels = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)
  for(r in 1:n_r){
    for(s in 1:n_s){
      fit = nlminb(start = init, objective = fit.index, lower = lower, upper = upper,
                   ind.obs = ind.vast[,s,r], N_3d = N_3d, F_3d = F_3d, M_m = M_m, survey_frac_v = survey_frac_v,
                   survey_sel_switch = survey_sel_switch, survey_sel_pars = survey_sel_pars,
                   dome_pars = dome_pars, size_bm = size_bm,  years_count = years_count, n_l = n_l, s=s, r=r, # this makes it so we can loop through season and region
                   control = list(eval.max = 500, iter.max = 500, rel.tol = 1e-8))
      q.vec[s,r] <- fit$par
      ind.pred.temp[,,s,r] <- N_3d[,,s,r]*exp(-survey_frac_v[s]*(F_3d[,,s]+M_m[,,s]))*q.vec[s,r]*survey.sels
      ind.pred[,s,r] <- apply(ind.pred.temp[,,s,r],1,sum)
    }
  }
  return(ind.pred)
}

############### generate data (with random error) ##########################
# Adds lognormal observation error
# Appropriate for catch data and indices of abundance
add.error.lognormal <- function(data,cv)
{
  #dim.data   = dim(data)  # row - # of year; col - # of data set (fleet or survey)
  #data.error = matrix(NA,nrow=dim.data[1],ncol=dim.data[2])

  #for (n in 1:dim.data[2])
  #{
  sd               = sqrt(log((cv^2+1)))
  data.error       = rlnorm(length(data),log(data),sd)
  #}
  return(data.error)
}

# For adding error to length composition information - not currently in use

# add.error.multino <- function(data,ess)
# {
#   dim.data      = dim(data)  # row - # of year; col - # of data set (fleet or survey)
#   data.error    = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
#   sample.temp   = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
#
#   original.com=sweep(data,1,rowSums(data),"/")
#
#   for (i in 1:dim.data[1])
#   {
#     if (sum(data[i,])>0)
#     {
#       sample.temp[i,] = rmultinom(1, ess[i], prob=original.com[i,])
#       data.error[i,]  = sample.temp[i,]/sum(sample.temp[i,])
#     }
#   }
#   return(data.error)
# }

# These functions used to generate the .dat files needed for the EM
# I don't think we need these anymore since we can just use the same
# data structures as the OM when we build projections and the MSE wrapper

# catch.data.gen <- function(year_start,year_end,catch.obj,catch.cv,catch.ess,add.error)
# {
#
#   n.season       = lapply(catch.obj, dim)$catch_4d[3]
#   n.fleets       = lapply(catch.obj, dim)$catch_4d[4]
#   n.lengthbin    = lapply(catch.obj, dim)$catch_4d[2]
#   n.year         = lapply(catch.obj, dim)$catch_4d[1]
#
#   n_s = n.season
#   years_count = n.year
#   n_fleets = n.fleets
#
#   year_start = year_start
#   year_end = year_end
#   catch.data.m   = matrix(-1,nrow=years_count*n_s*n_fleets,ncol=9+length(size_bm)) #changed n.lengthbins to what I've been using to count them 'size_bm'
#   catch.data.m
#   length(size_bm)
#
#   year      = rep(seq(year_start,year_end,1),n_s*n_fleets)
#   season    = rep(rep(1:n_s,each=years_count),n_fleets)
#   fleet     = rep(rep(1:n_fleets,each=years_count*n_s))
#   catch.tot = c()
#   cv        = c()
#   cpue      = rep(-1,years_count*n_s*n_fleets)
#   cpue.flag = rep(1,years_count*n_s*n_fleets)
#   cpue.cv   = rep(-1,years_count*n_s*n_fleets)
#   ess       = c()
#   comp      = data.frame()
#
#   for (f in 1:n_fleets)
#   {
#     for (s in 1:n_s)
#     {
#
#       cv.temp          = catch.cv[,f,s]
#       cv               = c(cv,cv.temp)
#       cv
#       catch.tot.temp   = catch.obj$catch_3d[,f,s]
#       catch.tot        = c(catch.tot,catch.tot.temp)
#       catch.tot
#       ess.temp         = catch.ess[,f,s]
#       ess              = c(ess,ess.temp)
#       ess
#       comp.temp        = catch.obj$catch_4d[,,s,f]
#       comp             = rbind(comp,comp.temp)
#       comp
#
#       if(add.error==TRUE)
#       {
#         catch.tot    = add.error.lognormal(catch.tot,cv)
#         comp         = add.error.multino(comp,ess)
#       }
#
#     }
#   }
#
#   catch.data.m[,1]                   = year
#   catch.data.m[,2]                   = season
#   catch.data.m[,3]                   = fleet
#   catch.data.m[,4]                   = catch.tot
#   catch.data.m[,5]                   = cv
#   catch.data.m[,6]                   = cpue
#   catch.data.m[,7]                   = cpue.flag
#   catch.data.m[,8]                   = cpue.cv
#   catch.data.m[,9]                   = ess
#   catch.data.m[,10:(9+n.lengthbin)]  = as.matrix(comp)
#   catch.data.m
#   return(catch.data.m)
#
# }


# # survey.data.gen <- function(survey_periods,survey_months,n_survey,survey.obj,survey.cv,survey.ess,add.error)
# # {
# #   nrow           = 0
# #   nrow.temp      = 0
# #   year           = c()
# #   indexN         = c()
# #   month          = c()
# #   index          = c()
# #   cv             = c()
# #   ess            = c()
# #   comp           = data.frame()
# #   n_sizeb        = lapply(survey.obj, dim)[[1]][2] #1 gives dimensions, 2 calls the columns' dimensions
# #   n_sizeb
# #   for (i in 1:n_survey)
# #   {
# #     nrow.temp    = lapply(survey.obj, dim)[[i]][1]
# #     nrow         = nrow + nrow.temp
# #     nrow
# #     year.temp    = seq(survey_periods[i,1],survey_periods[i,2],1)
# #     year         = c(year,year.temp)
# #     year
# #     length(year)
# #     indexN.temp  = rep(i,nrow.temp)
# #     indexN       = c(indexN,indexN.temp)
# #     indexN
# #     length(indexN)
# #     month.temp   = rep(survey_months[i],nrow.temp)
# #     month        = c(month,month.temp)
# #     month
# #     length(month)
# #     index.temp   = survey.obj[[i]][,1]
# #     index        = c(index,index.temp)
# #     index
# #     length(index)
# #     cv.temp      = survey.cv[[i]]
# #     cv           = c(cv,cv.temp)
# #     cv
# #     ess.temp     = survey.ess[,i,] #Hot fix gets it to work. I think the main error here has to do with the array being years_count rows
# #     ess          = c(ess,ess.temp)
# #     ess
# #     comp.temp    = survey.obj[[i]][,2:n_sizeb]
# #     comp         = rbind(comp,comp.temp)
# #     comp
# #     colSums(comp)
# #     dim(comp)
# #
# #     if(add.error==TRUE)
# #     {
# #       index      = add.error.lognormal(index,cv)
# #       comp       = add.error.multino(comp,ess)
# #
# #     }
# #
# #   }
# #
# #   survey.data.m  = matrix(-1,nrow=length(index),ncol=(6+n_sizeb-1))
# #   survey.data.m
# #   survey.data.m[,1]                 = year
# #   survey.data.m[,2]                 = indexN
# #   survey.data.m[,3]                 = month
# #   survey.data.m[,4]                 = index
# #   survey.data.m[,5]                 = cv
# #   survey.data.m[,6]                 = ess
# #   survey.data.m[,7:(6+n_sizeb-1)]   = as.matrix(comp)
# #   survey.data.m
# #   return(survey.data.m)
# #
# # }
#
# survey.data.gen <- function(survey_periods,survey_months,n_survey,survey.obj,survey.cv,survey.ess,add.error)
# {
#   nrow           = 0
#   nrow.temp      = 0
#   year           = c()
#   indexN         = c()
#   month          = c()
#   index          = c()
#   cv             = c()
#   ess            = c()
#   comp           = data.frame()
#   n_sizeb        = lapply(survey.obj, dim)[[1]][2]
#
#
#   for (i in 1:n_survey)
#   {
#     nrow.temp    = lapply(survey.obj, dim)[[i]][1]
#     nrow         = nrow + nrow.temp
#     year.temp    = seq(survey_periods[i,1],survey_periods[i,2],1)
#     year         = c(year,year.temp)
#     indexN.temp  = rep(i,nrow.temp)
#     indexN       = c(indexN,indexN.temp)
#     month.temp   = rep(survey_months[i],nrow.temp)
#     month        = c(month,month.temp)
#     index.temp   = survey.obj[[i]][,1]
#     index        = c(index,index.temp)
#     cv.temp      = survey.cv[[i]]
#     cv           = c(cv,cv.temp)
#     ess.temp     = survey.ess[,i,]
#     ess          = c(ess,ess.temp)
#     comp.temp    = survey.obj[[i]][,2:n_sizeb]
#     comp         = rbind(comp,comp.temp)
#
#     if(add.error==TRUE)
#     {
#       index      = add.error.lognormal(index,cv)
#       comp       = add.error.multino(comp,ess)
#     }
#
#   }
#   nrow
#   year
#   indexN
#   month
#   index
#   cv
#   ess
#   comp
#
#
#   survey.data.m  = matrix(-1,nrow=length(index),ncol=(6+n_sizeb-1))
#
#   survey.data.m[,1]                 = year
#   survey.data.m[,2]                 = indexN
#   survey.data.m[,3]                 = month
#   survey.data.m[,4]                 = index
#   survey.data.m[,5]                 = cv
#   survey.data.m[,6]                 = ess
#   survey.data.m[,7:(6+n_sizeb-1)]   = as.matrix(comp)
#
#   return(survey.data.m)
#
# }
#
