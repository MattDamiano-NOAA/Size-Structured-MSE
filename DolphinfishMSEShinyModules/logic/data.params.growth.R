########### Growth #################
# Von Bertalanffy growth function parameters
# Set up minimum and maximum size and increments

Linc = 100 # Increment by which size changes, in mm, # Has important consequences for the growth transition matrix
Lmin = size_bm[1]
Lmax = size_bm[20]

# values below based on values combined by sex from Schwenke and Buckel (2008)
Linf = 1290 # 7/17/2020: Does not have to be lower than Lmax
SELinf = 25.95 # Standard error of Linf # CV = sd/mean, sd = CV*mean
vbcv = 0.2 # CV for VBGF # vbcv = SELinf/Linf # pretty precise if we do this! we don't actually have this par # placeholder value, doesn't seem to matter - see below
K = 1.27 # K
SEK = 0.08 # Standard error of K
rhoLinfK = 0.5 # correlation between Linf and K, we don't know this, but we do know that the parameters are highly correlated
Linf_v = matrix(Linf,nrow = years_count, ncol = 1) # for time-varying growth, create matrices for VBGF parameters - probably won't use
SELinf_v = matrix(SELinf, nrow = years_count, ncol = 1)
K_v = matrix(K, nrow = years_count, ncol = 1)
SEK_v = matrix(SEK, nrow = years_count, ncol = 1)
rhoLinfK_v = matrix(rhoLinfK, nrow = years_count, ncol = 1)
n_growblock = n_block # set growth blocks to time-varying block
growblock_m <- matrix(NA, nrow = years_count, ncol = n_growblock) # matrix of growblocks for data generation
growth_array <- get_GM(size_bins_count,years_count,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
# check matrix
# growth_array[,,1]
# check that all columns sum to 1
# rowSums(growth_array[,,1])

# Something interesting - since we are getting a G for a year, fish in the first
# bin are most probable to grow a lot, and that makes biological sense
# Does this mean there's something off with the temporal structure or that growth is just that fast?
# it is sensitive to rholinfk, somewhat
# pretty insensitive to the CV of the VBGF
