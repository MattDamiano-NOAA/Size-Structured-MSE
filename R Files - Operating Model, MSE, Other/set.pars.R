# Size-structured MSE 
# AM parameter estimation file
# M. Damiano
# Last updated: 8-23-2023

# Set up priors/penalties for estimated quantities in the AM
# Set fourth value in each vector to 1 or -1 to turn "on" or "off"

###### Natural mortality priors ######
M_phase = -1 # set to "off"

M_prior = matrix(c(0.2, 0.1, 0.5, M_phase, 0.2, 0.0, 1), nrow = n_s, ncol = 7, byrow = TRUE) 

# Lorenzen natural mortality b parameter prior
M_b_prior = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = n_s, ncol = 7, byrow = TRUE)

# M_devs_priors
M_devs_phase = -1
M_devs_sd = 0.2
M_devs_lambda = 0

M_devs_prior = c(0, -10, 10, M_devs_phase, M_devs_sd, M_devs_lambda, 3)
M_devs_sd_prior = c(M_devs_sd, 0.0001, 10, -1, 1, M_devs_lambda, 1)

###### Fishing mortality priors ######
# cpue_power_prior
cpue_power_prior = c(1, 0.5, 1.5, -1, 0.5, 0, 1)
# selectivity prior
f_sel_pars_prior = matrix(c(10, 0, 1000, 1, 1, 0, 1), nrow = 8, ncol = 7, byrow = TRUE)

f_first_prior = matrix(c(1, 0, 10, 1, 1, 0.0, 1), nrow = n_s*n_fleets, ncol = 7, byrow = TRUE)

f_devs_prior = matrix(c(0, -1000, 1000, 1, -1, 0.0, 1), nrow = n_s*n_fleets, ncol = 7, byrow = TRUE)

###### Growth Priors ######
Linf_prior = c(500, 0, 800, -1, 0.5, 0, 1)
K_prior = c(0.173, 0, 0.5, -1, 0.5, 0, 1)
sd_Linf_prior = c(0.1, 0, 0.5, -1, 0.5, 0, 1)
sd_K_prior = c(0.1, 0, 0.5, -1, 0.5, 0, 1)
rho_prior = c(0.1, 0, 0.5, -1, 0.5, 0, 1)

###### Population dynamics priors ######
N_init_prior = c(N_init_tot, 0, 5e+10, 1, 1e+06, 0, 1)

###### Recruitment priors ######
# R proportion prior
R_l_pro_prior = c(0.1, 0, 0.5, -1, 0.5, 0, 1)
# R-S priors
alpha_prior = c(alpha, alpha*0.001, alpha*1000, -1, 0.2, 0, 1)
beta_prior = c(beta, beta*0.01, beta*100, -1, 0.2, 0, 1)
# environmental effects prior
R_env_prior = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = n_env, ncol = 7, byrow = TRUE)
# rec autocorrelation prior
R_autocor_prior = c(0, 0, 1, -1, -1, 0, 1)
# # rec devs prior information
R_devs_v_prior = c(0, -10, 10, 1, 1, 0, 1)
sd_R_devs_prior = c(sd_Rdevs, 0, 1, 1, -1, 0, 1)

###### Sex-change priors ######
Rsex = c(5, 0, 100, -1, 0.5, 0, 1)
# # length at 50% selection prior
L50_prior = c(200, 0, 500, -1, 0.5, 0, 1)

###### Survey priors ######
# # survey catchability power prior
q_power_prior = c(1, 0, 1.5, -1, 0.5, 0, 1)
# # survey selectivity prior
survey_sel_prior = matrix(c(20, 0, 1000, -1, 1, 0.0, 1), nrow = n_survey*2, ncol = 7, byrow = TRUE)

# Number of years projected
num.proj = 1
