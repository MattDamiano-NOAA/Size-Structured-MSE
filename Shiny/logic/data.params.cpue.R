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
