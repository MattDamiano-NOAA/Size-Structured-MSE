########### Maturity #################
# Define maturity parameters (for calculating SSB, not for yield)
L50 = 465 # based on avg of male and female L50 from Schwenke & Buckel (2008), Dolphinfish fork length @ 50% maturity
inter = 6.5 # same as ^

# set maturity parameters - need to re-name later
fe_prop_pars_m = matrix(NA, nrow = years_count, ncol = 2)
for (n in 1:years_count) {
  fe_prop_pars_m[n,] = c(L50, inter)
}
