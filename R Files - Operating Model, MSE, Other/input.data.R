# Size-structured MSE 
# Input data file
# M. Damiano
# Last updated: 8-23-2023

###### Model dimensions ######

terminal_year = 2021

# set the number of years 
years = rep(1990:terminal_year)

# set annual time step
n_t = length(years)

year_start = years[1]

year_end = tail(years, n = 1)

# set seasonal timestep; 1=no seasonality
n_s = c(1,2)[1] # can include however many seasons you want
# set blocks for time-varying quantities, parameters
# set to annual
n_block = 1

# growth increments 
size_bins = c(80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500)

# midpoint of growth increments
size_bm = c(90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490)

# set length structure variable based on number of size bins
n_l = length(size_bm)

###### Biological specifications #######

# set weight-length relationship 
WL_pars_m = matrix(NA, nrow = n_t, ncol = 2)
for (n in 1:n_t)
{
  WL_pars_m[n,]=c(-9.9,2.8) 
}

# Define maturity/sex change variables
L50 = 180 # length at 50% maturity
F50 = 120 # length at 50% female maturity (for sex changing species only)
inter = quantile(size_bm, 0.75)-quantile(size_bm,0.25) # interquartile range for use in the maturity function
# set maturity/sex change parameters
fe_prop_pars_m = matrix(NA, nrow = n_t, ncol = 2)
for(n in 1:n_t)
{
  fe_prop_pars_m[n,]=c(L50, inter)
}
# data.gen.R file will need matrices of these values to create the assessment input files
L50_mat = matrix(c(180), ncol = 1, nrow = n_t, byrow = TRUE)
maturity_m <- get_matu_prop(fe_prop_pars_m,F50,size_bm,n_t)

# Set up fecundity
# can be used with fecundity function if enabled in the Functions.R file
# fec_a = 7.69 # placeholder values for realistic relationship
# fec_b = 0.0053
# fec_pars_m = matrix(NA, nrow = n_t, ncol = 2)
# for(n in 1:n_t)
# {
#   fec_pars_m[n,]=c(fec_a, fec_b)
# }
# fec_mat <- get_fec(fec_pars_m, size_bm, weight_m, n_t)

####### Natural mortality #######
M_year_v = matrix(0.375, nrow = n_t, ncol = 1) 
# matrix for size effect on M
M_size_v = matrix(1, nrow = 1, ncol = n_l) 
# Fit to data for Mdevs (2), estimate (1) or turn off (0)
M_devs_flag = 0
# for parameters.dat file
M_value_as = 0.375
M_m <- get_M(M_size_v,M_year_v)

###### Selectivity #######
# Logistic selectivity parameters
sel_pars_v = c(300, 0.1) 
# Dome-shaped selectivity parameters
dome_pars = c(125,0.2,280,0.3) 

# set gear selectivity curve type
# 2 = Logistic
# 3 = Dome-shaped
sel_switch = 2
# Set up switch vector with selectivities for each fleet, e.g., 2, 2, 3 for logistic commercial and recreational, dome-shaped discard
switch_v = matrix(data=c(2,2,3), nrow = 3, ncol=1) 

####### Fishing Mortality ######

# Choose F values for aggregated fleets: commercial, recreational, discard
# F values are log-scale sums by fleet from 2023 black sea bass stock assessment

# Estimated F by fleet during 1990-2021 
# Commercial
# commercial hook and line
F.cl <- c(0.066, 0.066, 0.056, 0.046, 0.055, 0.040, 0.037, 0.047, 0.064, 0.059, 0.029, 0.027, 0.030, 0.025, 0.031, 0.023,
          0.022, 0.019, 0.020, 0.027, 0.022, 0.014, 0.026, 0.048, 0.073, 0.044, 0.053, 0.059, 0.051, 0.052, 0.031, 0.046)
F.cl <- as.numeric(F.cl)
# commercial pot
F.cp <- c(0.103, 0.090, 0.084, 0.078, 0.092, 0.076, 0.089, 0.096, 0.081, 0.121, 0.091, 0.111, 0.093, 0.096, 0.133, 0.096, 
          0.119, 0.087, 0.085, 0.121, 0.086, 0.069, 0.045, 0.052, 0.036, 0.041, 0.029, 0.068, 0.074, 0.082, 0.042, 0.025)
F.cp <- as.numeric(F.cp)
# avg commercial F
sum.comm.F <- array(NA, dim = c(31, 1))
for(t in 1:n_t){
  sum.comm.F[t] <- sum(F.cl[t], F.cp[t]) # sums all commercial F
}
# Recreational
# headboat 
F.hb <- c(0.045, 0.035, 0.027, 0.019, 0.020, 0.019, 0.021, 0.022, 0.021, 0.038, 0.025, 0.031, 0.021, 0.021, 0.041, 0.035, 
          0.032, 0.043, 0.026, 0.038, 0.066, 0.050, 0.023, 0.023, 0.021, 0.019, 0.019, 0.019, 0.028, 0.029, 0.027, 0.028)
F.hb <- as.numeric(F.hb)
# General private (MRIP)
F.rec <- c(0.107, 0.139, 0.135, 0.115, 0.165, 0.098, 0.172, 0.132, 0.088, 0.147, 0.208, 0.339, 0.185, 0.173, 0.446, 0.380, 
           0.295, 0.644, 0.555, 0.411, 0.794, 0.585, 0.388, 0.640, 1.411, 0.779, 0.611, 1.024, 0.541, 0.814, 0.633, 1.045)
F.rec <- as.numeric(F.rec)
sum.rec.F <- c()
for(t in 1:n_t){
  sum.rec.F[t] <- sum(F.hb[t], F.rec[t]) # sum all recreational F
}
# Discards
# comm discard
F.cd <- c(0.001, 0.001, 0.001, 0.001, 0.002, 0.001, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
          0.001, 0.002, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.001, 0.0, 0.0)
F.cd <- as.numeric(F.cd)
# hb discard
F.hbd <- c(0.0, 0.005, 0.001, 0.001, 0.003, 0.001, 0.002, 0.002, 0.001, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.0,
           0.001, 0.002, 0.002, 0.003, 0.003, 0.005, 0.006, 0.006, 0.007, 0.008, 0.010, 0.011, 0.012, 0.017, 0.014, 0.026)
F.hbd <- as.numeric(F.hbd)
# rec discard
F.recd <- c(0.009, 0.013, 0.013, 0.016, 0.026, 0.012, 0.015, 0.019, 0.015, 0.025, 0.035, 0.033, 0.024, 0.025, 0.047, 0.040,
            0.047, 0.067, 0.063, 0.051, 0.062, 0.080, 0.094, 0.075, 0.206, 0.178, 0.201, 0.322, 0.211, 0.337, 0.348, 0.391)
F.recd <- as.numeric(F.recd)
# avg disc F
sum.disc.F <- c()
for(t in 1:n_t){
  sum.disc.F[t] <- sum(F.cd[t], F.hbd[t], F.recd[t]) # can sum all discard F, but this ends up smoothing out the dome in this example
}
# initial F
F_init_v <- c(sum.comm.F[1],sum.rec.F[1],sum.disc.F[1])
n_fleets = length(F_init_v)
# create empty array for F by year, season and fleet
F_array = array(NA,c(n_t,n_s,n_fleets)) 
F_array[1:n_t,,1] <- sum.comm.F
F_array[1:n_t,,2] <- sum.rec.F
F_array[1:n_t,,3] <- sum.disc.F # need to set to rec discards to capture pattern
F_check1 <- F_array # check
# selectivity blocks and parameters for data generation files
f_sel_n_bloc = n_fleets
f_sel_switch_v = matrix(NA, ncol = 1, nrow = n_fleets)
for(i in 1:n_fleets)
{
  f_sel_switch_v[i] = switch_v[i]
}
f_sel_block_array = array(1, dim=c(n_t,n_fleets,n_s))
f_sel_block_array[,1,1] = rep(1,n_t)
f_sel_block_array[,2,1] = rep(2,n_t)
f_sel_block_array[,3,1] = rep(3,n_t)

###### Growth ######
# Von Bertalanffy growth function parameters
# Set up minimum and maximum size and increments 
Linc = 20
Lmin = size_bm[1]
Lmax = size_bm[21]
# Increment by which size changes
# Has important consquences for the growth transition matrix. 

# 7/17/2020: Does not have to be lower than Lmax
Linf = 510
# CV for VBGF
vbcv = 0.3
# Standard error of Linf
SELinf = vbcv*Linf
#K 
K = 0.183
# Standard error of K
SEK = vbcv*K
# correlation between Linf and K
rhoLinfK = 0.7
# for time-varying growth, create matrices for VBGF parameters
Linf_v = matrix(Linf,nrow = n_t, ncol = 1)
SELinf_v = matrix(SELinf, nrow = n_t, ncol = 1)
K_v = matrix(K, nrow = n_t, ncol = 1)
SEK_v = matrix(SEK, nrow=n_t, ncol = 1)
rhoLinfK_v = matrix(rhoLinfK, nrow = n_t, ncol=1)
# set growth blocks to time-varying block
n_growblock = n_block
# matrix of growblocks for data generation (not currently in use)
growblock_m <- matrix(NA, nrow = n_t, ncol = n_growblock)
growth_array <- get_GM(n_l,n_t,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)

####### Recruitment #######
# Define recruitment parameters
alpha = exp(18.076) # mean recruitment
beta = 400e+10 # stock size (eggs) at which we achieve half of max recruitment
# SR relationship, 2 = 2 parameter BH
indicator = 1 
# Spawning month
ssb_month = 3 # March for black sea bass

####### Population Dynamics #######

# initial population size
N_init_tot = 1.3e+08 
# proportions of N that go to each size bin
N_init_comp = c(0.3, 0.3, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005, 0.005, 0.005, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025) # changed 5/16/2022
# N_init_comp=N_init_comp/sum(N_init_comp) # if rescaling is required

# create vector for recruitment deviations
R_devs_v = matrix(NA, nrow = n_t, ncol = 1)

# log-scale recdevs from 2023 assessment during 1990-2021
assmt.rec.devs <- c(0.135034200463, -0.247454055582, -0.277203438790, 0.0905731492523, -0.00832834312386, -0.111905195994, -0.0548709812026, 0.272250224791, 0.00105824077504,
                    0.218993963583, 0.159496743529, 0.131538958968, 0.0677757483408, 0.246053687813, -0.0461199762185, 0.133100977517, 0.155505026508, 0.232091366233, 
                    0.311432460601, 0.622424128633, 0.347560525606, 0.0149544395932, -0.0492959345655, -0.150851488311, -0.557937781735, -0.902317223341, -0.906884180543,
                    -1.10110882787, -1.53287911936, -1.85179985638, -1.142154, -1.142154)

# set proportion of recruits going to each length class, depends on how we define recruitment
R_l_pro = c(0.5,0.25,0.15,0.10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
# number of length classes to which fish recruit
n_R_l = length(which(R_l_pro != 0))
# sd for rec devs
sd_Rdevs = 0.536 
# scale recruitment deviations if necessary
for(t in 1:n_t){
  R_devs_v[t] <- assmt.rec.devs[t]
}
R_devs_v <- scale(R_devs_v, center = T, scale = F) 

# fit environmental effects to rec devs? 1 = yes, 0 = no
R_fit_env = 0 # not currently in use, but environmental effects can be linked to recruitment this way

# number of environmental data sets
n_env = 3

# environmental data matrix (year x n_env)
env_data = matrix(0, nrow = n_t, ncol = n_env)

###### Survey/Index ######

# set number of surveys
n_survey = 1 
# set mean month when survey occurs; this can be a vector containing months for each survey
survey_months = 5
# set survey catchability 
survey_q = 1E-6

# select survey selectivity parameters
survey_sel_pars = c(160, 0.3)
# matrix of survey selectivity parameters
survey_sel_pars_m = matrix(NA, ncol = 2, nrow = n_t)
for (n in 1:n_t)
{
  survey_sel_pars_m[n,] = c(160, 0.3) 
}

# select survey selectivity relationship; 2 = logistic
survey_sel_switch = 2 
# generate survey period (time series of when survey took place)
survey_periods = matrix(NA, ncol = 2, nrow = n_survey)
for(n in 1:n_survey)
{
  survey_periods[n,] = c(1990, terminal_year)
}

###### Set up observation error #######

# effective sample size for length composition data
ess = 100 
# cv for catch
cv = 0.05
# cv for survey
s.cv = 0.27
# create arrays of catch cv and ess by year, fleet and season
catch.cv = array(cv, c(n_t, n_fleets, n_s))
catch.ess = array(ess, c(n_t, n_fleets, n_s))
# create arrays of survey cv and ess by year, survey and season
survey.cv = array(s.cv, c(n_t, n_s, n_survey))
survey.ess = array(ess, c(n_t, n_s, n_survey))

