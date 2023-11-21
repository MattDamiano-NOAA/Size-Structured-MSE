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
survey_sels <- cal_sels(survey_sel_pars,survey_sel_switch,dome_params,size_bm)

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
ind.vast <- array(NA,dim = c(years_count,seasons_count,regions_count))
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
catch.cv = array(cv, c(years_count, n_fleets, seasons_count))
catch.ess = array(ess, c(years_count, n_fleets, seasons_count))
# generate data with error: set to TRUE or FALSE
# add.error = FALSE # This needs to be commented out during MSE loop (MD 12-23-21)
# create arrays of survey cv and ess by year, survey and season
survey.cv = array(s.cv, c(years_count, seasons_count, n_survey))
survey.ess = array(ess, c(years_count, seasons_count, n_survey))
