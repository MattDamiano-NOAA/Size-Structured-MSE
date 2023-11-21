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
