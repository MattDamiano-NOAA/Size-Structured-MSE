########### Natural Mortality #################
# Molto et al. 2022 suggest M can be as high as 0.25 per month (1.00 per quarter?)
# set M to seasonal value of 1.00 for now
# annual M of 4.00 isn't far off from what Oxenford (1999) report, which is around 3.5
M_year_v = matrix(1, nrow = years_count, ncol = 1) # keep these dimensions and just apply M each season - same effect
M_size_v = matrix(1, nrow = 1, ncol = size_bins_count) # matrix for size effect on M
