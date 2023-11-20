# Dolphinfish MSE
# Dolphinfish OM: Input data file
# Last modified: 10-11-2023 by M. Damiano

initial_year <- 1986
terminal_year <- 2022

years = rep(initial_year:terminal_year) # set the number of years
years_count = length(years) # set annual time step
year_start = years[1] # set start year
year_end = tail(years, n = 1) # a bit redundant, but there it is

n_s = 4 # set seasonal timestep; 1=no seasonality
n_block = 1 # set blocks for time-varying quantities, parameters; 1=no blocks

n_r = 7 # number of regions for spatially-explicit model (we may want 8 to split the Caribbean into north and south later)
size_bins = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000) # define size bins
size_bm = c(150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550, 1650, 1750, 1850, 1950, 2050) # midpoints
n_l = length(size_bm) # set length structure variable based on number of size bins

# Biological specifications

########### Weight-Length #################
# set weight-length relationship parameters
WL_pars_m = matrix(NA, nrow = years_count, ncol = 2)
for (n in 1:years_count)
{
  WL_pars_m[n,] = c(log(2e-8),2.8) # parameters an avg based on several studies outlined in Oxenford (1999)
}

########### Maturity #################
# Define maturity parameters (for calculating SSB, not for yield)
L50 = 465 # based on avg of male and female L50 from Schwenke & Buckel (2008), Dolphinfish fork length @ 50% maturity
inter = 6.5 # same as ^

# set maturity parameters - need to re-name later
fe_prop_pars_m = matrix(NA, nrow = years_count, ncol = 2)
for (n in 1:years_count) {
  fe_prop_pars_m[n,] = c(L50, inter)
}


########### Natural Mortality #################
# Molto et al. 2022 suggest M can be as high as 0.25 per month (1.00 per quarter?)
# set M to seasonal value of 1.00 for now
# annual M of 4.00 isn't far off from what Oxenford (1999) report, which is around 3.5
M_year_v = matrix(1, nrow = years_count, ncol = 1) # keep these dimensions and just apply M each season - same effect
M_size_v = matrix(1, nrow = 1, ncol = n_l) # matrix for size effect on M

########### Selectivity #################
# choose selectivity parameters

# fleets with logistic selectivity - based on length frequency data
# selectivity = likelihood of catching certain indivs during a fishing operation
log_pars = c(100, 0.1) # dummy when not using logistic pars
pll_comm_sel = c(800, 0.005) # based on US PLL (comm) - U.S. Pelagic Long Line commercial selectivity
priv_rec_sel_S = c(740, 0.01) # based on PR & FL length dist (MRIP)

# fleets with dome-shaped selectivity
dome_pars = c(400,0.2,1200,0.3) # dummy when not using dome-shaped pars
priv_rec_sel_N = c(300, 0.01, 1100, 0.01) # patterns are very similar for N states priv and overall for hire
# fh_rec_sel = c(400, 0.01, 1200, 0.01)
disc_gen_sel = c(100,0.01,500,0.01) # general discard selectivity

# 10-10-2023: 6 fleets for now: commercial, private rec N & S, for hire N & S, and general discard
# No spatial structure for now; can use areas as fleets approach later
sels_4d = array(NA,c(years_count,n_l,n_s,6)) # haven't defined 'n_fleets' yet
for (s in 1:n_s) {
  for (t in 1:years_count) {
    sels_4d[t,,s,1] = cal_sels(pll_comm_sel, 2, dome_pars, size_bm) # commercial
    sels_4d[t,,s,2] = cal_sels(log_pars, 3, priv_rec_sel_N, size_bm) # private N
    sels_4d[t,,s,3] = cal_sels(priv_rec_sel_S, 2, dome_pars, size_bm) # private S
    sels_4d[t,,s,4] = cal_sels(log_pars, 3, priv_rec_sel_N, size_bm) # for-hire N
    sels_4d[t,,s,5] = cal_sels(log_pars, 3, priv_rec_sel_N, size_bm) # for-hire S (same pattern as N)
    sels_4d[t,,s,6] = cal_sels(log_pars, 3, disc_gen_sel, size_bm) # general dead discard fleet
  }
}

########### Fishing Mortality #################
# Choose F values for fleets
# bear in mind that these are summed to single values

# For dolphinfish, we'll need an F for each fleet, but not necessarily by area
# We can use selectivity to modify F based on area (areas as fleets approach)

# create empty array for F by year, season and fleet
# dummy values for now for one comm, one rec, one discard
# What fisheries do we know we need coastwide?
# Commercial (maybe combined or just longline is ok);
# Private recreational x 2
# For-hire recreational x 2
# High seas intl (eventually)
# Discard fleet
# Then we will probably need a Venezuelan artisanal drift gill net for lower Caribbean
F_init_v <- c(0.1, 0.3, 0.3, 0.05, 0.05, 0.05)
# Fleets: commercial, private rec N, private rec S, for-hire N, for-hire S, discard
n_fleets = length(F_init_v)

F_array = array(NA,c(years_count,n_s,n_fleets))
F_array[1:years_count,,1] <- F_init_v[1]
F_array[1:years_count,,2] <- F_init_v[2]
F_array[1:years_count,,3] <- F_init_v[3]
F_array[1:years_count,,4] <- F_init_v[4]
F_array[1:years_count,,5] <- F_init_v[5]
F_array[1:years_count,,6] <- F_init_v[6]
# F_check1 <- F_array

########### Growth #################
# Von Bertalanffy growth function parameters
# Set up minimum and maximum size and increments

Linc = 100
Lmin = size_bm[1]
Lmax = size_bm[20]
# Increment by which size changes
# Has important consequences for the growth transition matrix

# 7/17/2020: Does not have to be lower than Lmax
# values below based on values combined by sex from Schwenke and Buckel (2008)
Linf = 1290
# Standard error of Linf
# CV = sd/mean, sd = CV*mean
SELinf = 25.95
# CV for VBGF
# vbcv = SELinf/Linf # pretty precise if we do this! we don't actually have this par
vbcv = 0.2 # placeholder value, doesn't seem to matter - see below
#K
K = 1.27
# Standard error of K
SEK = 0.08
# correlation between Linf and K, we don't know this, but we do know that the parameters are highly correlated
rhoLinfK = 0.5
# for time-varying growth, create matrices for VBGF parameters - probably won't use
Linf_v = matrix(Linf,nrow = years_count, ncol = 1)
SELinf_v = matrix(SELinf, nrow = years_count, ncol = 1)
K_v = matrix(K, nrow = years_count, ncol = 1)
SEK_v = matrix(SEK, nrow = years_count, ncol = 1)
rhoLinfK_v = matrix(rhoLinfK, nrow = years_count, ncol = 1)
# set growth blocks to time-varying block
n_growblock = n_block
# matrix of growblocks for data generation
growblock_m <- matrix(NA, nrow = years_count, ncol = n_growblock)
growth_array <- get_GM(n_l,years_count,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
# check matrix
growth_array[,,1]
# check that all columns sum to 1
rowSums(growth_array[,,1])

# Something interesting - since we are getting a G for a year, fish in the first
# bin are most probable to grow a lot, and that makes biological sense
# Does this mean there's something off with the temporal structure or that growth is just that fast?
# it is sensitive to rholinfk, somewhat
# pretty insensitive to the CV of the VBGF

########### Recruitment #################
# Define recruitment parameters
alpha = 10000000 # we don't actually know this number
R_devs_v = matrix(NA, nrow = years_count, ncol = 1) # deviations in mean recruitment
sd_recdevs = 1
rec.sto = TRUE # stochastic mean recruitment?
if (rec.sto == FALSE){
  alpha = alpha
  for (t in 1:years_count){
    R_devs_v[t] <- 0
  }
} else {
  ln.rec <- rnorm(1, mean = log(alpha), sd = sd_recdevs)
  alpha = exp(ln.rec)
  for (t in 1:years_count){
    R_devs_v[t] <- rnorm(1, 0, 1)
  }
}

# beta = 400e+10 #stock size (eggs) at which we achieve half of max recruitment; not used in mean recruitment model
# SR relationship, 2 = 2 parameter BH
indicator = 1 # Beverton-Holt
# Spawning month: if they spawn all year, then this might not be needed
ssb_month = 3
# Call SSB fractions
# SSB_frac_v <- get_ssb.fracs(ssb_month, n_s)

########### Population Dynamics #################
# Create initial abundance proportions at size
# 6/25/2020: The first size class, the one to which recruits recruit, should
# be zero in the initial population by size vector.
N_init_tot = 3000000 # slightly higher than max estimated catch in private rec fishery MRIP time series

# initial composition is arbitrary for now, will probably need multiple hypotheses to quantify uncertainty
#100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
N_init_comp = c(0.2, 0.125, 0.1, 0.1, 0.1, 0.1, 0.050, 0.050, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.01, 0.01, 0.01, 0.01, 0.01)
# check that it sums to 1.0
sum(N_init_comp)

# set proportion of recruits going to each length class, depends on how we define recruitment
R_l_pro = c(0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0,0,0,0,0,0,0)
# check length
# length(R_l_pro)
# # check that it sums to 1
# sum(R_l_pro)
# number of length classes to which fish recruit
n_R_l = length(which(R_l_pro != 0)) # not sure this is used anymore

# R_devs_v <- scale(R_devs_v, center = T, scale = F)

# EM remnant, but we may still try to fit rec devs to environmental indices, e.g., tripole (Volkov et al. 2019)
# for environmental effects
# fit environmental effects to rec devs? 1 = yes, 0 = no
R_fit_env = 0

# number of environmental data sets
n_env = 3

# environmental data matrix (year x n_env)
# consider putting tripole data here for use later
env_data = matrix(0, nrow = years_count, ncol = n_env)

########### Effort (for CPUE calcs) #################
# test case is private rec effort, which looks logarithmic
# Still need to think about how to deal with seasonality
# these trends are annual - for pll, they are the same regardless of season (surprisingly)

# Probably the best approach is to make assumptions about effort like this,
# and then look at it without effort information to see if it even matters (to the MSE and MP)

# pars in each simulated effort dataset come from fits in Excel (i.e., tested add trendline until I got the best r-squared)
x <- seq(1,years_count,1)

# polynomial  private rec north of FL (effort in directed trips)
p_rec_eff_N <- -1251.1*(x^2) + 48606*x + 1e06
# plot(p_rec_eff_N, ylim = c(0,2500000))

# polynomial private rec in FL and east (includes PR)
p_rec_eff_S <- -22.6*(x^2) + 5253.3*x - 18633
# plot(p_rec_eff_S, ylim = c(0,250000)) # check

# logarithmic for-hire north of FL (effort in directed trips)
p_x <- log(x)
fh_rec_eff_N <- 3692.2*p_x - 1432.3 # this is the simplest thing to do: take the slope and int from Excel fits
# plot(fh_rec_eff_N, ylim = c(0,20000)) # check

# polynomial for-hire rec in FL and south (includes PR)
fh_rec_eff_S <- 112.04*(x^2) + -5355.3*x + 104038
# plot(fh_rec_eff_S, ylim = c(0,120000)) # check

# polynomial for commercial pll (effort in 1000s of hooks)
pll_comm_eff <- -2.74*(x^2) + 48*x + 4675
# plot(pll_comm_eff, ylim = c(0,9000)) # check

# assume discards track with N of FL rec for now

# Need a switch to apply different parameterizations of functions to get CPUE
# match with F_init_v
# commercial, private rec N, private rec S, for-hire N, for-hire S, discard
# 1 = polynomial for commercial (pll)
# 2 = polynomial for private rec N
# 3 = polynomial for private rec S
# 4 = logarithmic for for-hire N
# 5 = polynomial for for-hire S
# assume discard follows private rec S where smaller fish are not legal to keep
#
eff_switch_v <- c(1,2,3,4,5,3) # set it to match up with F_init_v

# make an array to input to the function
eff_array <- array(NA, dim = c(years_count,n_fleets))
eff_array[,1] <- pll_comm_eff
eff_array[,2] <- p_rec_eff_N
eff_array[,3] <- p_rec_eff_S
eff_array[,4] <- fh_rec_eff_N
eff_array[,5] <- fh_rec_eff_S
eff_array[,6] <- p_rec_eff_S

########### Survey #################
# set number of surveys
n_survey = 4 # 4 "surveys" per year for each region
# set mean month when survey occurs; this can be a vector containing months for each survey
survey_months_v = c(2,5,8,11) # "survey" or VAST index for now, is reported quarterly
# set survey catchability
#survey_q = exp(-7.67)
survey_q_v = rep(1E-6,7) # starting q is the same for all regional indices for now

# select survey selectivity parameters
survey_sel_pars = c(500, 0.3) # assuming knife-edged logistic sel with L50 = 500 mm, approx 20 inches
# matrix of survey sel pars
survey_sel_pars_m = matrix(NA, ncol = 2, nrow = years_count)
for (n in 1:years_count)
{
  survey_sel_pars_m[n,] = c(500, 0.3)
}
# Note: it looks like the survey selectivity is 1, which might make sense?
# select survey selectivity relationship; 2 = logistic
survey_sel_switch = 2
# generate survey period (time series of when survey took place)
survey_periods = matrix(NA, ncol = 2, nrow = n_survey)
for (n in 1:n_survey)
{
  survey_periods[n,] = c(initial_year, terminal_year)
}

# call survey fractions
survey_frac_v <- c()
survey_frac_v[1] = (survey_months_v[1] + 1)/12
survey_frac_v[2] = (survey_months_v[2] + 1)/12
survey_frac_v[3] = (survey_months_v[3] + 1)/12
survey_frac_v[4] = (survey_months_v[4] + 1)/12
survey_sels <- cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)

# read in VAST indices
d <- read.csv('data/VASTindex.csv')

# d <- as.data.frame(cbind(ind.vast$NED,ind.vast$VBMN,ind.vast$NNC,ind.vast$NCFL,ind.vast$SAR,ind.vast$FLK,ind.vast$CAR))
# This isn't enough - you need to get it into a 37 x 4 x 7 array
# That means sorting by season
library(dplyr)
library(tidyr)
ind.vast <- as_tibble(d)
winter <- d %>% filter(Season == "Winter") # apparently there's an extra space in the .csv
winter[is.na(winter)] <- 0
spring <- d %>% filter(Season == "Spring")
spring[is.na(spring)] <- 0
summer <- d %>% filter(Season == "Summer")
summer[is.na(summer)] <- 0
fall <- d %>% filter(Season == "Fall")
fall[is.na(fall)] <- 0

# decided to brute force it because more like tedi(ous)verse, amirite?
ind.vast <- array(NA,dim = c(years_count,n_s,n_r))
# Set northern high seas index
ind.vast[,1,1] <- fall$NED
ind.vast[,2,1] <- winter$NED
ind.vast[,3,1] <- spring$NED
ind.vast[,4,1] <- summer$NED
# Set VA to Montauk index
ind.vast[,1,2] <- fall$VBMN
ind.vast[,2,2] <- winter$VBMN
ind.vast[,3,2] <- spring$VBMN
ind.vast[,4,2] <- summer$VBMN
# Set northern NC index
ind.vast[,1,3] <- fall$NNC
ind.vast[,2,3] <- winter$NNC
ind.vast[,3,3] <- spring$NNC
ind.vast[,4,3] <- summer$NNC
# Set southern NC to E. FL index
ind.vast[,1,4] <- fall$NCFL
ind.vast[,2,4] <- winter$NCFL
ind.vast[,3,4] <- spring$NCFL
ind.vast[,4,4] <- summer$NCFL
# Set Sargasso Sea/high seas index
ind.vast[,1,5] <- fall$SAR
ind.vast[,2,5] <- winter$SAR
ind.vast[,3,5] <- spring$SAR
ind.vast[,4,5] <- summer$SAR
# Set FL Keys index
ind.vast[,1,6] <- fall$FLK
ind.vast[,2,6] <- winter$FLK
ind.vast[,3,6] <- spring$FLK
ind.vast[,4,6] <- summer$FLK
# Set Caribbean Sea index
ind.vast[,1,7] <- fall$CAR
ind.vast[,2,7] <- winter$CAR
ind.vast[,3,7] <- spring$CAR
ind.vast[,4,7] <- summer$CAR

#Check and rescale
ind.vast <- ind.vast

# for testing q search bounds
init <- 1e-08
upper <- 1.0
lower <- 1e-12

########### Error #################
# set effective sample size
ess = 100
# set catch cv
cv = 0.05
# set index cv
s.cv = 0.27
# set catch time series starting and ending year
# year_start = 1978
# year_end = 2017
# create arrays of catch cv and ess by year, fleet and season
catch.cv = array(cv, c(years_count, n_fleets, n_s))
catch.ess = array(ess, c(years_count, n_fleets, n_s))
# generate data with error: set to TRUE or FALSE
# add.error = FALSE # This needs to be commented out during MSE loop (MD 12-23-21)
# create arrays of survey cv and ess by year, survey and season
survey.cv = array(s.cv, c(years_count, n_s, n_survey))
survey.ess = array(ess, c(years_count, n_s, n_survey))


