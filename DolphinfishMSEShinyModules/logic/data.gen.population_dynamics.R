# Dolphinfish MSE
# Dolphinfish OM: data generation file
# Last modified: 9-19-2023 by M. Damiano

weight_matrix    <- get_weight(WL_pars_m,size_bm,years_count)                    # Call weight-length relationship function to get weight matrix: weight_matrix
maturity_matrix  <- get_matu_prop(fe_prop_pars_m,size_bm,years_count)            # Call maturity function to get maturity matrix: maturity_m
M_matrix         <- get_M(M_size_v,M_year_v)                             # Call natural mortality function to get M matrix: M_m
f_object         <- get_F(F_array,sels_4d,n_fleets,size_bm,n_s,years_count,n_l)  # Call F functions and check F_4d and F_3d objects, F_inits/deviations
growth_array     <- get_GM(n_l,years_count,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
SSB_frac_v       <- get_ssb.fracs(ssb_month, n_s) # growth_array
PopDy            <- collect_population_dynamics(
  N_init_tot, N_init_comp,
  f_object$F_3d, M_matrix, SSB_frac_v, maturity_matrix,
  weight_matrix, alpha, beta, indicator,
  R_devs_v, R_l_pro, years_count, n_s, n_r, n_growblock,
  growth_array, mov.mat, mov.mat2, mov.mat3, mov.mat4
)
