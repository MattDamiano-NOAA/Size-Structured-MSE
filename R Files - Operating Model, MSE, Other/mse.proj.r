# Size-structured MSE 
# MSE projection function file
# M. Damiano
# Last updated: 8-23-2023

mse.proj <- function(n_proj,N_terminal,R_proj,sels_4D,F_apex_proj,M_vec_proj,G_proj,F_allocation_proj){ # probably needs to be a function of sels_4D
  
  N_proj <- array(NA, c(n_proj,n_l,n_s))  # projected abundance, no seasonality
  sels_proj <- array(NA, c(n_proj, n_l, n_s, n_fleets))
  F_proj_array <- array(NA, c(n_proj, n_l, n_s, n_fleets))
  M_matrix = matrix(NA, nrow=n_proj,ncol=n_l)
  
  SSB_proj <- c(n_proj) # projected SSB
  
  for (f in 1:n_fleets)
  {
    for (t in 1:n_proj)
    {
      for (s in 1:n_s)
      {
        sels_proj[t,,s,f] = sels_4D[t,,s,f]
        F_proj_array[t,,s,f]=F_allocation_proj[t,s,f]*sels_proj[t,,s,f] # dimensions of sels need to change for time varying selectivity here and in AM
      }
    }
  }
  
  F_matrix <- rowSums(F_proj_array, dims = 2)
  
  for (t in 1:n_proj){ # total F is used here because selectivity is the same among fleets 
    M_matrix[t,] = M_vec_proj
    if(t==1){
      N_proj[t,,1] = N_terminal%*%GM+R_proj[t]*R_l_pro
    }else{
      N_proj[t,,1] = (N_proj[t-1,,1]*exp(-(F_matrix[t,]+M_vec_proj)))%*%G_proj+R_proj[t]*R_l_pro
    }
    SSB_proj[t] = sum(N_proj[t,,1]*exp(-SSB_frac_v[1]*(F_matrix[t,]+M_vec_proj))*weight_m[1,]*maturity_m[1,])
  }
  
  C_proj_4d <- array(NA, c(n_proj, n_l, n_s, n_fleets)) # Projected catch
  C_proj_3d <- array(NA, c(n_proj, n_fleets, n_s))
 
  for (f in 1:n_fleets)
  {
    for (t in 1:n_proj)
    {
      for (s in 1:n_s)
      {
        Z = F_matrix[t,]+M_vec_proj
        C_proj_4d[t,,s,f] = N_proj[t,,1]*(1-exp(-Z))*((sels_proj[t,,s,f]*F_allocation_proj[t,s,f])/Z)
        C_proj_3d[t,f,s] = sum(C_proj_4d[t,,s,f])
      }
    }
  }
  # end catch
  
  # Get survey(s)
  survey.temp  <- array(NA,c(n_proj,n_l,n_survey))
  survey <- array(NA, c(n_proj, n_l+1, n_survey))
  
  for (n in 1:n_survey)
  {
    for (t in 1:n_proj)
    {
      if(n_s <= 1)
      {
        Z = F_matrix[t,]+M_vec_proj
        survey.sels.proj = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)
        survey.temp[t,,n] = N_proj[t,,1]*exp(-survey_frac_m[n,1]*Z)*survey_q[n]*survey.sels.proj
      }
      if(n_s > 1)
      {
        #temp = survey_frac_m[n,2]
        #survey.temp[t,,n] = N_proj[t,,temp]*exp(-survey_frac_m[n,3]*(F_sim + M_sim))*survey_q[n]*survey_sels_est[n,]
      }
    }
    survey[,,n] = cbind(apply(survey.temp[,,n],1,sum),survey.temp[,,n]) 
  }
  
  projections <- list(N_proj=N_proj, SSB_proj=SSB_proj, R_proj=R_proj, Catch4d = C_proj_4d, Catch3d = C_proj_3d, Survey = survey, F_proj = F_matrix, M_proj = M_matrix)                
  return(projections)
}
