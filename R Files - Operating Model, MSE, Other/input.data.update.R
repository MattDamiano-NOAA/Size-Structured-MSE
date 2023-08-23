# Size-structured MSE 
# Input data update file
# M. Damiano
# Last updated: 8-23-2023

library(abind) # used to update F_array object

###### set up weight-length parameters for projected years ######
# set weight-length relationship parameters
WL_pars_m = matrix(NA, nrow = n_t, ncol = 2)
for (n in 1:n_t)
{
  WL_pars_m[n,]=c(-9.9,2.8)
}
# maturity parameters for projected years
fe_prop_pars_m = matrix(NA, nrow = n_t, ncol = 2)
for(n in 1:n_t)
{
  fe_prop_pars_m[n,]=c(L50, inter)
}

# Fecundity parameters for projected years
# fec_pars_m = matrix(NA, nrow = n_t, ncol = 2)
# for(n in 1:n_t)
# {
#   fec_pars_m[n,]=c(7.69, 0.0053)
# }

# still need these for the data.gen.proj.R file
L50_mat = matrix(c(180), ncol = 1, nrow = n_t, byrow = TRUE) 
maturity_m <- get_matu_prop(fe_prop_pars_m,F50,size_bm,n_t)
M_year_v = matrix(0.375, nrow = n_t, ncol = 1) 

# Update F_array
F_proj.array <- array(NA, c(n_proj, n_s, n_fleets))
F_proj.v <- array(NA, c(1,n_fleets))
for(f in 1:n_fleets){
  F_proj.v[f] <- F_apex_proj*F_fl_prop[f]
}

for(f in 1:n_fleets){
  F_proj.array[,,f] <- F_proj.v[f]
}
F_array <- abind(F_array, F_proj.array, along = 1) 

# selectivity blocks and parameters for data generation files
f_sel_n_bloc = n_fleets
#f_sel_pars_m = sel_pars_m
f_sel_switch_v = matrix(NA, ncol = 1, nrow = n_fleets)
for(i in 1:n_fleets)
{
  f_sel_switch_v[i] = switch_v[i]
}
f_sel_block_array = array(1, dim=c(n_t,n_fleets,n_s))
f_sel_block_array[,1,1] = rep(1,n_t)
f_sel_block_array[,2,1] = rep(2,n_t)
f_sel_block_array[,3,1] = rep(3,n_t)

f_obj <- get_F(F_array,sel_pars_v,sel_switch,n_fleets,size_bm,n_s,n_t,n_l)
F_3d <- f_obj$F_3d
F_4d <- f_obj$F_4d
f_inits.obj <- get_F_inits(F_array)

# for time-varying growth, create matrices for VBGF parameters
Linf_v = matrix(Linf,nrow = n_t, ncol = 1)
SELinf_v = matrix(SELinf, nrow = n_t, ncol = 1)
K_v = matrix(K, nrow = n_t, ncol = 1)
SEK_v = matrix(SEK, nrow=n_t, ncol = 1)
rhoLinfK_v = matrix(rhoLinfK, nrow = n_t, ncol=1)
# set growth blocks to time-varying block
n_growblock = n_block
# matrix of growblocks for data generation
# not totally certain what this needs to contain...
growblock_m <- matrix(NA, nrow = n_t, ncol = n_growblock)
growth_array <- get_GM(n_l,n_t,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)

# Updated alpha term for alternate recruitment scenario projections
proj_alpha = mean(PopDy$Recruits[25:30])
alpha = proj_alpha
# alpha = alpha

R_devs_v.proj = matrix(NA, nrow = n_proj, ncol = 1)
for (t in 1:n_proj)
{
  # R_devs_v.proj[t] <- log(Stats.proj$R_proj[t])-log(alpha) # zero deviation from mean, no stochasticity
  R_devs_v.proj[t] <- rnorm(1, mean = 0, sd = sd_Rdevs) #sd based on sd of last 10 10 years
}
# turning this off in case it's contributing to non-convergence 10-11-2022
# R_devs_v.proj <- R_devs_v.proj - mean(R_devs_v.proj) # centers random lognormal deviates on 0
R_devs_v = rbind(R_devs_v,R_devs_v.proj)
R_devs_v = scale(R_devs_v)
env_data = matrix(0, nrow = n_t, ncol = n_env)

# ########### Survey #################
# # set number of surveys
# n_survey = 1 
# # set mean month when survey occurs; this can be a vector containing months for each survey
# survey_months = 6
# # set survey catchability 
# survey_q = exp(-7.67)
# # select survey selectivity parameters
# survey_sel_pars = c(3.41, 1.93)
# # matrix of survey sel pars
# survey_sel_pars_m = matrix(NA, ncol = 2, nrow = n_t)
# for (n in 1:n_t)
# {
#   survey_sel_pars_m[n,] = c(3.41, 1.93) # 6/25/2020: reversed these parameters for the same reason: I think AM reads them this way
# }
# # Note: it looks like the survey selectivity is 1, which might make sense? 
# # select survey selectivity relationship; 2 = logistic
# survey_sel_switch = 2 
# # generate survey period (time series of when survey took place)
# survey_periods = matrix(NA, ncol = 2, nrow = n_survey)
# for(n in 1:n_survey)
# {
#   survey_periods[n,] = c(1990, terminal_year)
# }
# 
# # call survey fractions
# survey_frac_m <- get_survey.fracs(survey_months,n_survey,n_s)
# 
# ########### Error #################
# # set effective sample size
# ess = 100 
# # set catch cv
# cv = 0.27 
# # set catch time series starting and ending year
# # year_start = 1978
# # year_end = 2017
# # create arrays of catch cv and ess by year, fleet and season
# catch.cv = array(cv, c(n_t, n_fleets, n_s))
# catch.ess = array(ess, c(n_t, n_fleets, n_s))
# # generate data with error: set to TRUE or FALSE
# add.error = FALSE
# # create arrays of survey cv and ess by year, survey and season
# survey.cv = array(cv, c(n_t, n_s, n_survey))
# survey.ess = array(ess, c(n_t, n_s, n_survey))










