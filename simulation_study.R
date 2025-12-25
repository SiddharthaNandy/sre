sre_parallel = function( N ){
  
  # Arguments -
  # N is the parallel distribution ID.  
  
  # Returns -
  # the N^th ID results.
  
  # rm( list = ls() )
  
  require(spatstat) # To use 'rjitter' in 'spatial_locations.R' and 'true_knots.R'
  require(sp) # To call 'spJitter' from pedometrics
  require(pedometrics) # To use 'spJitter'
  require(matrixcalc)
  require(splus2R) # To use 'vecnorm' 
  require(Rcpp) # To call '.cpp' functions
  require(MASS)
  require(Matrix)
  require(pracma)
  require(psych)
  require(spam) # For sparse matrices
  require(cIRT) # To use 'riwishart' in 'gen_data.R'

  # Sourcing all related functions ( Change the directory to the folder where you save all the .R files.) :
  # -------------------------------------------------------------------------------------------------------
  
  folder = "/mnt/ffs24/home/nandysid/SRE_selection/"
  
  # Construction of Spatial locations, MR knots, and basis function -
  source(paste( folder , "spatial_locations.R" , sep = "") )
  source(paste( folder , "mr_knots.R", sep = "") )
  source(paste( folder , "mr_radial_basis.R", sep = "") )
  source(paste( folder , "mr_basis.R", sep = "") )

  # Data generation for the SRE model -
  source(paste( folder , "true_knots.R", sep = "") )
  source(paste( folder , "gen_data.R", sep = "") )
  source(paste( folder , "cholesky_embedding.R", sep = "") )
  
  # Likelihood - Gradient - Fisher Information - Armijo line search - Group LASSO optimization -
  sourceCpp(paste( folder , "norm_cholesky.cpp", sep = "") )
  sourceCpp(paste( folder , "Sigma_inv_det.cpp", sep = "") )
  sourceCpp(paste( folder , "like_grad_info_j.cpp", sep = "") )
  sourceCpp(paste( folder , "BCGD_projection_j.cpp", sep = "") )
  sourceCpp(paste( folder , "armijo_line_update.cpp", sep = "") )
  source(paste( folder , "sre_grplasso.R", sep = "") )
  # sourceCpp(paste( folder , "sre_grplasso.cpp", sep = "") )
  
  # SRE Estimation -
  # source(paste( folder , "truep_falsep.R", sep = "") )
  source(paste( folder , "sre_estimation.R", sep = "") )
  
  # Values to test the function with -
  # N = 98
  
  T1 = Sys.time()
  
  D0 = 2
  n0 = 400 # choices from ( 100 , 225 , 400 , 625 , 900 )
  Rep0 = 75 # choices from ( 1 , 25 , 75 )
  Irreg0 = 0.009 # Fixed
  B_type = "W" # choices from ("LB","W")
  nres0 = 2 # number of resolutions
  NC0 = 3 # number of knots in one direction of the coarsest resolution
  NCb0 = 0 # number of buffer knots in one direction of the coarsest resolution
  knots = 5*NC0^2 - 4*NC0 + 1 # number of knots from two resolutions
  TrueM = 1 # (1 or 2)
  rn0 = 18 # (for Kn = 13, use rn - 4 , 6 , 8 , and Kn = 34, use rn = 14, 16, 18) true rank
  Iter0 = 50 # Vary this to improve more accuracy
  
  # Spatial locations :
  # -------------------
  
  sp_coords = spatial_locations( n = n0 , Irreg = Irreg0 )
  
  SRE_Armijo = sre_estimation( MC = N , D = D0 , n = n0 , Rep = Rep0 , coords = sp_coords , 
                               basis_type = B_type , nres = nres0 , NC = NC0 , NCb = NCb0 , 
                               Model = TrueM , rank = rn0 , Iter = Iter0 )
  
  fname = paste( folder , "n" , n0 , "T" , Rep0 , "Kn" , knots , "rn" , rn0 , "Model" , TrueM , B_type , "MC" , N , '.RD' , sep = '' ) 
  save( SRE_Armijo , file = fname ) 
  
  T2 = Sys.time()
  Ctime = T2 - T1  
  print( Ctime )
  print( c(N,n0,Rep0,rn0) )
  
  # result = paste(as.character(Sys.info()["nodename"]),' sre_parallel:',as.character(N), "\n")
  
}