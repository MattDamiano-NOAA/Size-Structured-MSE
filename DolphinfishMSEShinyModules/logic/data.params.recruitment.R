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
# SSB_frac_v <- get_ssb.fracs(ssb_month, seasons_count)
