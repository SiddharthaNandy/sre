spatial random effects selection

Simulation Study :
Users should call simulation_study.R to run simulation replication.

List of functions used in simulation_study.R -

  # Construction of Spatial locations, MR knots, and basis function -
  spatial_locations.R; 
  mr_knots.R; 
  mr_radial_basis.R; 
  mr_basis.R
  
  # Data generation for the SRE model -
  true_knots.R;
  gen_data.R; 
  cholesky_embedding.R; 
    
  # Likelihood - Gradient - Fisher Information - Armijo line search - Group LASSO optimization -
  norm_cholesky.cpp; 
  Sigma_inv_det.cpp; 
  like_grad_info_j.cpp; 
  BCGD_projection_j.cpp; 
  armijo_line_update.cpp; 
  sre_grplasso.R; 
  
  # SRE Estimation -
  sre_estimation.R

Real Data study:
Users should call sre_RD_run.R to run simulation replication.

List of functions used in sre_RD_run.R -

  # Construction of MR knots, and basis function -
  mr_knots.R
  mr_radial_basis.R
  mr_basis.R
      
  # Likelihood - Gradient - Fisher Information - Armijo line search - Group LASSO optimization -
  norm_cholesky.cpp
  Sigma_inv_det.cpp
  like_grad_info_j.cpp
  BCGD_projection_j.cpp
  armijo_line_update.cpp
  sre_grplasso_data.R
  
  # SRE Estimation -
  sre_DOI_WR.R


