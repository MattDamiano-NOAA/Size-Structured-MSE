# Dolphinfish MSE
# Dolphinfish OM: data generation file
# Last modified: 9-19-2023 by M. Damiano

weight_m    <- get_weight(WL_pars_m,size_bm,n_t) # Call weight-length relationship function to get weight matrix: weight_m
maturity_m  <- get_matu_prop(fe_prop_pars_m,size_bm,n_t) # Call maturity function to get maturity matrix: maturity_m
M_m         <- get_M(M_size_v,M_year_v) # Call natural mortality function to get M matrix: M_m
f_obj       <- get_F(F_array,sels_4d,n_fleets,size_bm,n_s,n_t,n_l) # Call F functions and check F_4d and F_3d objects, F_inits/deviations
F_3d        <- f_obj$F_3d
F_4d        <- f_obj$F_4d
f_inits.obj <- get_F_inits(F_array)
GM          <- cal_GM(Lmin,Lmax,Linc,Linf,SELinf,K,SEK,rhoLinfK) # Call growth transition matrix function and check growth_array
grow_3d     <- get_GM(n_l,n_t,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
growth_array = grow_3d # which name it has doesn't really matter, as long as it's used when calling the PopDy function
SSB_frac_v  <- get_ssb.fracs(ssb_month, n_s) # growth_array
PopDy       <- pop_dyn(N_init_tot,N_init_comp,F_3d,M_m,SSB_frac_v,maturity_m,weight_m,alpha,beta,indicator,R_devs_v,R_l_pro,n_t,n_s,n_r,n_growblock,growth_array,mov.mat,mov.mat2,mov.mat3,mov.mat4)
N_3d        <- PopDy$Abundance # PopDy is the object containing recruitment, abundance, and SSB
catch.obj   <- cal_catch(f_obj,M_m,N_3d,n_t,n_l,n_s,n_fleets,n_r) # Check catch object
catch_4d    <- catch.obj$catch_4d # catch.obj contains catch at size by fleet
catch_3d    <- catch.obj$catch_3d # and total catch (size-aggregated) by fleet
cpue_array  <- get_cpue(catch_3d, eff_switch_v, eff_array)
index       <- get.fitted.pred.index(ind.vast,N_3d,F_3d,M_m,survey_frac_v,survey_q_v,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,n_t,n_l,n_r)
