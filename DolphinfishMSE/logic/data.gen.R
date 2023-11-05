# Dolphinfish MSE
# Dolphinfish OM: data generation file
# Last modified: 9-19-2023 by M. Damiano

# Call OM functions to produce values for data gen
# Call weight-length relationship function to get weight matrix: weight_m
weight_m <- get_weight(WL_pars_m,size_bm,n_t)
# plot(weight_m[1,])
# Call maturity function to get maturity matrix: maturity_m
maturity_m <- get_matu_prop(fe_prop_pars_m,size_bm,n_t)
# plot(maturity_m[1,])
# Call natural mortality function to get M matrix: M_m
M_m <- get_M(M_size_v,M_year_v)
# Call F functions and check F_4d and F_3d objects, F_inits/deviations
f_obj <- get_F(F_array,sels_4d,n_fleets,size_bm,n_s,n_t,n_l)
F_3d <- f_obj$F_3d
F_4d <- f_obj$F_4d
f_inits.obj <- get_F_inits(F_array)
# Call growth transition matrix function and check growth_array
GM <- cal_GM(Lmin,Lmax,Linc,Linf,SELinf,K,SEK,rhoLinfK)
grow_3d <- get_GM(n_l,n_t,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
growth_array = grow_3d # which name it has doesn't really matter, as long as it's used when calling the PopDy function
# growth_array
SSB_frac_v <- get_ssb.fracs(ssb_month, n_s)
PopDy <- pop_dyn(N_init_tot,N_init_comp,F_3d,M_m,SSB_frac_v,maturity_m,weight_m,alpha,beta,indicator,R_devs_v,R_l_pro,n_t,n_s,n_r,n_growblock,growth_array,mov.mat,mov.mat2,mov.mat3,mov.mat4)
N_3d <- PopDy$Abundance # PopDy is the object containing recruitment, abundance, and SSB
# Check catch object
catch.obj <- cal_catch(f_obj,M_m,N_3d,n_t,n_l,n_s,n_fleets,n_r)
# catch.obj contains catch at size by fleet
catch_4d <- catch.obj$catch_4d
# and total catch (size-aggregated) by fleet
catch_3d <- catch.obj$catch_3d
# survey_frac_v <- get_survey.fracs(survey_months_v,n_survey,n_s)
# survey.obj <- cal_survey_inds(N_3d,F_3d,M_m,survey_frac_v,survey_q_v,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,n_t,n_l,n_r)
# survey.obj

# test obs error function
# data <- catch_3d
# cv = 0.3
# data.error <- add.error.lognormal(data,cv)
# # works great
# rel.err <- (data.error-catch_3d)/catch_3d
# mean(rel.err)

#
cpue_array <- get_cpue(catch_3d, eff_switch_v, eff_array)

index <- get.fitted.pred.index(ind.vast,N_3d,F_3d,M_m,survey_frac_v,survey_q_v,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,n_t,n_l,n_r)

# All of this below is old code for writing .dat files to be read into AD Model Builder software
# IGNORE

###################################################

# set file path
#input.dir <- "/Users/jcao22/Google Drive/Students/Damiano/MARFIN MSE/Assessment model and figure code/InputFiles"
#input.season.dir <- "/Users/jcao22/Google Drive/Students/Damiano/MARFIN MSE/Assessment model and figure code/InputFiles/Season"
#input.year.dir <- "/Users/jcao22/Google Drive/Students/Damiano/MARFIN MSE/Assessment model and figure code/InputFiles/Year"

# write control file
# sink(paste(input.dir,"/Control.DAT",sep=''))
# cat("#Indicator of time-step","\n")
# cat(c(n_s),"\n")
# cat("#Number of years","\n")
# cat(c(n_t),"\n")
# cat("#Number of seasons in each year","\n")
# cat(c(n_s),"\n")
# cat("#Number of months in each season","\n")
# cat(c(rep(12/n_s,n_s)),"\n")
# cat("#First year of the input data","\n")
# cat(c(year_start),"\n")
# cat("#First year for calculation","\n")
# cat(c(year_start),"\n")
# cat("#Last year for calculation","\n")
# cat(c(year_end_retro),"\n")
# cat("#Include likelihood constants? (1=yes; 0=no)","\n")
# cat(c(0),"\n")
# cat("#DiCohortTrackingBeginYear (track a cohort) ","\n")
# cat(c(year_start),"\n")
# cat("# Check Number","\n")
# cat(c(999),"\n")
# sink()
#
# # write Biology_Data.DAT
#
# sink(paste(input.dir,"/Biology_Data.DAT",sep=''))
# cat("#Number of size bins","\n")
# cat(n_l,"\n")
# cat("#Lower and upper boundary for each size bin (mm) ","\n")
# cat(size_bins,"\n")
# cat("#Parameters of weight-length relation by year (W in gram and L in mm) log(w)=a+b*log(L)//year a b","\n")
# for (i in 1:n_t){
#   cat(as.matrix(cbind(seq(year_start,year_end,1),WL_pars_m))[i,],"\n")
# }
# cat("#Maturity","\n")
# for (i in 1:n_t){
#   cat(as.matrix(cbind(seq(year_start,year_end,1),maturity_m))[i,],"\n")
# }
# # cat("#Natural mortality dev flag (0-turn off; 1-estimate; 2-fit to data)","\n")
# # cat(M_devs_flag,"\n")
# cat("#Natural mortality weighting factors by size ","\n")
# cat(M_size_v,"\n")
# cat("#Natural mortality weighting factors by year","\n")
# cat(rep(1,n_t),"\n")
# cat("#Number of size bins to which recruitment recruits","\n")
# cat(n_R_l,"\n")
# cat("#Proportions of the annual recruitment recruits to each season (Will not be used when the time-step is year)","\n")
# cat(1,"\n")
# cat("#Spawning month (the beginning of the month)","\n")
# cat(ssb_month,"\n")
# cat("#Recruitment-SSB relation","\n")
# cat(indicator,"\n")
# cat("#Env fit to Recruitment devs","\n")
# cat(R_fit_env,"\n")
# cat("#Number of environmental variables for RS (max=3;has to be 3)","\n")
# cat(n_env,"\n")
# cat("#Values of environmental variables (year*variable)","\n")
# for (i in 1:n_t){
#   cat(env_data[i,],"\n")
# }
# cat("#Initial condition set-up (8 choices: 0-7; see technical document)","\n")
# cat(1,"\n")
# cat("#Proportions-at-size (Pia vector; will only be used when initial condition is set to 1 or 2)","\n")
# cat(N_init_comp,"\n")
# cat("# Check Number","\n")
# cat(c(999),"\n")
# sink()
#
# # write catch.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/CatchDataYear.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/CatchDataSeason.DAT",sep=''))
#
# cat("# Number of fleets ","\n")
# cat(n_fleets,"\n")
# cat("# Catch unit for each fleet (0-Number(Million); 1-Biomass(1000mt))","\n")
# for (i in 1:n_s){
#   cat(rep(0,n_fleets),"\n")
# }
# cat("# Start size bin of selectivity for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(1,n_fleets),"\n")
# }
# cat("# End size bin of selectivity for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(n_l,n_fleets),"\n")
# }
# cat("# Likelihood function of length composition data for each fleet (1-multinomial; 2-robust normal)","\n")
# for (i in 1:n_s){
#   cat(rep(1,n_fleets),"\n")
# }
# cat("# Likelihood function of total catch for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(4,n_fleets),"\n")
# }
# cat("# Likelihood function of CPUE for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(4,n_fleets),"\n")
# }
# cat("# Lambda value of composition component in objective function for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(1,n_fleets),"\n")
# }
# cat("# Lambda value of total catch in objective function for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(1,n_fleets),"\n")
# }
# cat("# Lambda value of CPUE in objective function for each fleet","\n")
# for (i in 1:n_s){
#   cat(rep(0,n_fleets),"\n")
# }
# cat("# Number of data points for catch data ","\n")
# cat(nrow(catch.data),"\n")
# cat("#Year, timestep, Fleet, Total Catch ","\n")
# for (i in 1:nrow(catch.data)){
#   cat(catch.data[i,], "\n")
# }
# cat("# Use CPUE? (0-no; 1-yes)","\n")
# for (i in 1:n_s){
#   cat(rep(0,n_fleets),"\n")
# }
# cat("# Number of CPUE catchability  (Number of time blocks)","\n")
# cat(1,"\n")
# cat("# Option for CPUE catchability calculation method","\n")
# cat(1,"\n")
# cat("# CPUE catchability time blocks set-up","\n")
# for (s in 1:n_s){
#   for (i in 1:n_t){
#     cat(c(seq(year_start,year_end,1)[i],rep(1,n_fleets)),"\n")
#   }
# }
# cat("# Number of fleet selectivity time blocks","\n")
# cat(f_sel_n_bloc,"\n")
# cat("# Fleet selectivity option for each fleet","\n")
# cat(f_sel_switch_v,"\n")
# cat("# Fleet selectivity time blocks set-up","\n")
# for (s in 1:n_s){
#   for (i in 1:n_t){
#     cat(c(seq(year_start,year_end,1)[i],f_sel_block_array[i,,s]),"\n")
#   }
# }
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
# # write Survey_Data.DAT
#
# sink(paste(input.dir,"/Survey_Data.DAT",sep=''))
# cat("#Number of available survey indices","\n")
# cat(n_survey,"\n")
# cat("#Unit of each survey index (1-biomass; 0-numbers)","\n")
# cat(rep(0,n_survey),"\n")
# cat("#Start size bin of selectivity for each survey","\n")
# cat(rep(1,n_survey),"\n")
# cat("#End size bin of selectivity for each survey","\n")
# cat(rep(n_l,n_survey),"\n")
# cat("#Use Index? (1-yes; 0-no)","\n")
# cat(rep(1,n_survey),"\n")
# cat("#Likelihood function of length composition data for each survey (1-multinomial; 2-robust normal)","\n")
# cat(rep(1,n_survey),"\n")
# cat("#Likelihood function of index for each survey","\n")
# cat(rep(4,n_survey),"\n")
# cat("#Lambda value of composition component in objective function for each survey","\n")
# cat(rep(1,n_survey),"\n")
# cat("#Lambda value of index in objective function for each survey","\n")
# cat(rep(1,n_survey),"\n")
# cat("#Number of data points for survey indices","\n")
# cat(c(nrow(survey.data)),"\n")
# cat("#Year, IndexNum,IndexMonth(the end of the month), Index value, CV, ESS, Proportions at size","\n")
# for (i in 1:nrow(survey.data)){
#   cat(survey.data[i,],"\n")
# }
# cat("#Lambda value of sex change component in objective function","\n")
# cat(c(1),"\n")
# cat("#IndexNum for proportions of female at size","\n")
# cat(c(1),"\n")
# cat("#Proportions of female at size","\n")
# for (i in 1:n_t){
#   cat(maturity_m[i,],"\n")
# }
# cat("#Number of survey catchability (allow time-varying catchability)","\n")
# cat(c(n_survey),"\n")
# cat("#Option for catchability calculation method  (1-qI=sum(ln(Iobs/B^e1))/NUears; 2-qI=ln(sum(Iobs/B^e1))/NUears)","\n")
# cat(rep(2,n_survey),"\n")
# cat("#Index catchability time blocks set-up","\n")
# for (i in 1:n_t){
#   cat(c(seq(year_start,year_end,1))[i],seq(1,n_survey,1),"\n")
# }
# cat("#Use fleet selectivity? (negative value-not use; fleet number-use that particular fleet)","\n")
# cat(rep(-2,n_survey),"\n")
# cat("#Number of survey selectivity time blocks","\n")
# cat(n_survey,"\n")
# cat("#Survey selectivity option for each survey (1-by size; 2-logisitic; 3-double logistic; 4-double normal)","\n")
# cat(rep(2,n_survey),"\n")
# cat("#Survey selectivity time blocks set-up","\n")
# for (i in 1:n_t){
#   cat(c(seq(year_start, year_end,1))[i],seq(1,n_survey,1),"\n")
# }
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
# # write BPR.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/BRP_Data_Year.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/BRP_Data_Season.DAT",sep=''))
#
# cat("#Maximum F value in penalty term","\n")
# cat(5,"\n")
# cat("#Selectivity set-up (-1-input; 0-averaged fleet selectivity; #?#:fleet selectivity)","\n")
# cat(0,"\n")
# cat("#Selectivity input (if above option is set to -1)","\n")
# for (i in 1:n_s){
#   cat(sels,"\n")
# }
# cat("#Equilibrium period","\n")
# cat(50,"\n")
# cat("#Reference year for mortality","\n")
# cat(n_t,"\n")
# cat("#F proportions for each season (1 for annual time-step)","\n")
# cat(rep(1,n_s),"\n")
# cat("#Growth matrix set-up (which time block of growth matrix)","\n")
# cat(rep(1,n_s),"\n")
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
#
# # write GrowthMatrix.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/GrowthMatrix_Year.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/GrowthMatrix_Season.DAT",sep=''))
#
# cat("#Growth matrix set-up (2-one growth matrix for a year using growth fraction for each season;  1-VBGF Parameters; 0-input growth matrix)","\n")
# cat(1,"\n")
# cat("#Number of growth matrices","\n")
# cat(1,"\n")
# cat("#Growth matrix time blocks set-up","\n")
# for (i in 1:n_t){
#   cat(c(seq(1,n_t,1)[i],rep(1,1),"\n"))
# }
# cat("#Growth proportion for each Season","\n")
# cat(rep(n_s/12,n_s),"\n")
# # cat("#Growth Transition Matrix", "\n")
# # cat(GM, "\n")
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
#
# # write Parameters_Ini.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/Parameters_Ini_Year.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/Parameters_Ini_Season.DAT",sep=''))
#
# cat("#Fleet Selectivity Parameters","\n")
# # for (i in 1:f_sel_n_bloc){
# #   cat(f_sel_pars_m[i,],"\n") # needs to be changed such that fleet 3 has dome pars
# # }
# cat("#Fleet 1 selectivity pars","\n")
# cat(sel_pars_v, "\n")
# cat("#Fleet 2 selectivity pars","\n")
# cat(sel_pars_v, "\n")
# cat("#Fleet 3 selectivity pars","\n")
# # cat(sel_pars_v, "\n") # change to "dome pars" if you want dome-shaped
# cat(dome_pars, "\n")
# cat("#Fishing mortality of the first year for each fleet","\n")
# for (i in 1:n_fleets){
#   cat(f_inits.obj$log_f_first_year[i,],"\n")
# }
# cat("#Fishing mortality deviations for each year and fleet (fleet outer loop, year inner loop)","\n")
# cat("#Fleet 1", "\n")
# for (t in 1:n_t){
#     cat(f_inits.obj$log_F_devs[t,,1],"\n")
# }
# cat("#Fleet 2", "\n")
# for (t in 1:n_t){
#   cat(f_inits.obj$log_F_devs[t,,2],"\n")
# }
# cat("#Fleet 3", "\n")
# for (t in 1:n_t){
#   cat(f_inits.obj$log_F_devs[t,,3],"\n")
# }
# cat("#CPUE catchability power parameter for each time block","\n")
# cat(1,"\n")
# cat("#Survey Index Selectivity parameters","\n")
# for (i in 1:n_survey){
#   cat(survey_sel_pars_m[i,],"\n")
# }
# cat("#Survey catchability power parameter for each time block","\n")
# cat(rep(1,n_survey),"\n")
# cat("#Initial condition parameters","\n")
# cat(N_init_tot,"\n")
# cat("#R-S relationship parameters (2)","\n")
# cat(c(log(alpha),beta),"\n")
# cat("#R-S environment coefficients","\n")
# cat(rep(0,n_env),"\n")
# cat("#Recruitment deviation (log scale)","\n")
# cat(R_devs_v,"\n")
# cat("#Recruitment autocorrelation coefficient","\n")
# cat(0,"\n")
# cat("#Standard deviation of recruitment deviation in log scale","\n")
# cat(sd_Rdevs,"\n")
# cat("#Natural mortality","\n")
# cat(M_value_as,"\n")
# cat("#Lorenzen natural mortality b","\n")
# cat(rep(0,n_s),"\n")
# cat("#VBGF Parameters for each time block (Linf)","\n")
# cat(Linf,"\n")
# cat("#VBGF Parameters for each time block (K)","\n")
# cat(K,"\n")
# cat("#VBGF Parameters for each time block (SD of Linf)","\n")
# cat(SELinf,"\n")
# cat("#VBGF Parameters for each time block (SD of K)","\n")
# cat(SEK,"\n")
# cat("#VBGF Parameters for each time block  abs(Rho between Linf and K)","\n")
# cat(rhoLinfK,"\n")
# cat("#Proportion of recruitment-at-size","\n")
# cat(R_l_pro,"\n")
# cat("#L50 for each year","\n")
# for (i in 1:n_t){
#   cat(L50_mat[i],"\n")
# }
# cat("#Rsex","\n")
# cat(log(inter),"\n")
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
# # write Prior.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/Prior_Year.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/Prior_Season.DAT",sep=''))
#
# cat("#Natural mortality ","\n")
# for (i in 1:n_s){
#   cat(M_prior[i,],"\n")
# }
#
# cat("#Lorenzen natural mortality b","\n")
# for (i in 1:n_s){
#   cat(M_b_prior[i,],"\n")
# }
# # cat("#Natural M devs","\n")
# # cat(M_devs_prior,"\n")
# # cat("#Natural M_sd devs","\n")
# # cat(M_devs_sd_prior,"\n")
#
# cat("#R-S Relationship Parameters (2)","\n")
# cat(alpha_prior,"\n")
# cat(beta_prior,"\n")
# cat("#R-S Relationship Environment coefficients","\n")
# for (i in 1:n_env){
#   cat(R_env_prior[i,],"\n")
# }
# cat("#Recruitment autocorrelation coefficient","\n")
# cat(R_autocor_prior,"\n")
# cat("#Recruitment deviations (log scale)","\n")
# cat(R_devs_v_prior,"\n")
# cat("#Standard deviation of recruitment deviation in log scale","\n")
# cat(sd_R_devs_prior,"\n")
# cat("#CPUE catchability power parameter","\n")
# cat(cpue_power_prior,"\n")
# cat("#Fleet selectivity parameters","\n")
# for (i in 1:8){ # must match number of sel pars for all fleets
#   cat(f_sel_pars_prior[i,],"\n")
# }
# cat("#Survey catchability power parameter","\n")
# for (i in 1:n_survey){
#   cat(q_power_prior,"\n")
# }
# cat("#Survey selectivity parameters","\n")
# for (i in 1:(n_survey*2)){
#   cat(survey_sel_prior[i,],"\n")
# }
# cat("#Initial condition Parameters","\n")
# cat(N_init_prior,"\n")
# cat("#Fishing mortality of the first year for each fleet","\n")
# for (i in 1:(n_s*n_fleets)){
#   cat(f_first_prior[i,],"\n")
# }
# cat("#Fishing mortality deviation for each fleet","\n")
# for (i in 1:(n_s*n_fleets)){
#   cat(f_devs_prior[i,],"\n")
# }
# cat("#VBGF Parameters for each time block","\n")
# for (i in 1:n_growblock){
#   cat(Linf_prior,"\n")
#   cat(K_prior,"\n")
#   cat(sd_Linf_prior,"\n")
#   cat(sd_K_prior,"\n")
#   cat(rho_prior,"\n")
# }
# cat("#Proportion of recruitment-at-size","\n")
# cat(R_l_pro_prior,"\n")
# cat("#L50","\n")
# cat(L50_prior,"\n")
# cat("#Rsex","\n")
# cat(Rsex,"\n")
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
# # write Projection.DAT
#
# if (n_s==1)
#   sink(paste(input.year.dir,"/Projection_Year.DAT",sep=''))
# if (n_s>1)
#   sink(paste(input.season.dir,"/Projection_Season.DAT",sep=''))
#
# cat("#Do projections? (1=yes, 0=no)","\n")
# cat(0,"\n")
# cat("#Number of years for projection","\n")
# cat(num.proj,"\n")
# cat("#Recruitment set-up (1-Input recruitment; 2-Input recruitment deviation and use R-S relationship; 3-Randomly choose from history)","\n")
# cat(1,"\n")
# cat("#Recruitment input (if above option is set to 1, still need to enter values even if not 1)","\n")
# for(i in 1:num.proj){
#   cat(1000, "\n")
# }
# cat("#Selectivity option (1-Input selectivity; 2-BRP Selectivity; 3-Averaged values of fleet selectivity)","\n")
# cat(1,"\n")
# cat("#Selectivity input (if above option is set to 1, still need to enter values even if not 1)","\n")
# for(i in 1:n_s){
#   for(j in 1:1){
#     cat(rep(1,n_l),"\n")
#   }
# }
# cat("#Control? (1-Catch; 0-F)","\n")
# cat(0,"\n")
# cat("#Catch or F input","\n")
# for(i in 1:num.proj){
#   cat(rep(0.5,n_s),"\n")
# }
# cat("#Growth matrix set-up (which time block of growth matrix)","\n")
# for(i in 1:num.proj){
#   cat(rep(1,n_s),"\n")
# }
# cat("# Check Number","\n")
# cat(c(-22122),"\n")
# sink()
#
#
