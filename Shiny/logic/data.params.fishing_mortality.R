########### Fishing Mortality #################
# Choose F values for fleets
# bear in mind that these are summed to single values

# For dolphinfish, we'll need an F for each fleet, but not necessarily by area
# We can use selectivity to modify F based on area (areas as fleets approach)

# create empty array for F by year, season and fleet
# dummy values for now for one comm, one rec, one discard

F_array = array(NA,c(years_count,seasons_count,n_fleets))
F_array[1:years_count,,1] <- F_init_v[1]
F_array[1:years_count,,2] <- F_init_v[2]
F_array[1:years_count,,3] <- F_init_v[3]
F_array[1:years_count,,4] <- F_init_v[4]
F_array[1:years_count,,5] <- F_init_v[5]
F_array[1:years_count,,6] <- F_init_v[6]
