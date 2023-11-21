# Dolphinfish MSE
# Dolphinfish OM: data generation file
# Last modified: 9-19-2023 by M. Damiano

weight_matrix    <- get_weight(WL_pars_m,size_bm,years_count)                    # Call weight-length relationship function to get weight matrix: weight_matrix
maturity_matrix  <- get_matu_prop(fe_prop_pars_m,size_bm,years_count)            # Call maturity function to get maturity matrix: maturity_m
Mortality_matrix <- get_M(M_size_v,M_year_v)                                     # Call natural mortality function to get M matrix: Mortality_matrix
f_object         <- get_F(F_array,sels_4d,n_fleets,size_bm,seasons_count,years_count,size_bins_count)  # Call F functions and check F_4d and F_3d objects, F_inits/deviations
growth_array     <- get_GM(size_bins_count,years_count,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
SSB_frac_v       <- get_ssb.fracs(ssb_month, seasons_count) # growth_array
PopDy            <- collect_population_dynamics(
  N_init_tot, N_init_comp,
  f_object$F_3d, Mortality_matrix, SSB_frac_v, maturity_matrix,
  weight_matrix, alpha, beta, indicator,
  R_devs_v, R_l_pro, years_count, seasons_count, regions_count, n_growblock,
  growth_array, movement.matrix1, movement.matrix2, movement.matrix3, movement.matrix4
)
