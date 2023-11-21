########### Selectivity #################
# choose selectivity parameters

# fleets with logistic selectivity - based on length frequency data
# selectivity = likelihood of catching certain indivs during a fishing operation
logistic_params = c(100, 0.1) # dummy when not using logistic pars
selectivity.pelagic_long_line_commercial = c(800, 0.005) # based on US PLL (comm) - U.S. Pelagic Long Line commercial selectivity
selectivity.private_recreational_south = c(740, 0.01) # based on PR & FL length dist (MRIP)

# fleets with dome-shaped selectivity
dome_params = c(400,0.2,1200,0.3) # dummy when not using dome-shaped pars
selectivity.private_recreational_north = c(300, 0.01, 1100, 0.01) # patterns are very similar for N states priv and overall for hire
# selectivity.for_hire = c(400, 0.01, 1200, 0.01)
selectivity.general_discards = c(100,0.01,500,0.01) # general discard selectivity

# 10-10-2023: 6 fleets for now: commercial, private rec N & S, for hire N & S, and general discard
# No spatial structure for now; can use areas as fleets approach later
sels_4d = array(NA,c(years_count,size_bins_count,seasons_count,6)) # haven't defined 'n_fleets' yet
for (s in 1:seasons_count) {
  for (t in 1:years_count) {
    sels_4d[t,,s,1] = cal_sels(selectivity.pelagic_long_line_commercial, 2, dome_params, size_bm) # commercial
    sels_4d[t,,s,2] = cal_sels(logistic_params, 3, selectivity.private_recreational_north, size_bm) # private N
    sels_4d[t,,s,3] = cal_sels(selectivity.private_recreational_south, 2, dome_params, size_bm) # private S
    sels_4d[t,,s,4] = cal_sels(logistic_params, 3, selectivity.private_recreational_north, size_bm) # for-hire N
    sels_4d[t,,s,5] = cal_sels(logistic_params, 3, selectivity.private_recreational_north, size_bm) # for-hire S (same pattern as N)
    sels_4d[t,,s,6] = cal_sels(logistic_params, 3, selectivity.general_discards, size_bm) # general dead discard fleet
  }
}
