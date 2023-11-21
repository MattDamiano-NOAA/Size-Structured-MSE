
########### Weight-Length #################
# set weight-length relationship parameters
WL_pars_m = matrix(NA, nrow = years_count, ncol = 2)
for (n in 1:years_count)
{
  WL_pars_m[n,] = c(log(2e-8),2.8) # parameters an avg based on several studies outlined in Oxenford (1999)
}
