rm( list = ls() )
# require(parallel)

folder = "/home/sxn581/SRE_selection/"
setwd(folder)

# Sourcing the function :
# -----------------------

source("sre_DOI_WR.R")

load("coords_loc.RData")
load("DOI.RData")


require(splines) # To use 'bs'
require(fields) # To use the temperature dataset of Colorado
require(matrixcalc)
require(splus2R)
require(Rcpp) # To call '.cpp' files
require(MASS)
require(Matrix)
require(pracma)
require(psych)
require(spam) # For sparse matrices
require(cIRT) # To use 'riwishart' for starting value
require(autoFRK) # To compare with the other method with TPS basis
require(LatticeKrig) # To use MR knot location from the package and compare

# Sourcing all related functions ( Change the directory to the folder where you save all the .R files.) :
# -------------------------------------------------------------------------------------------------------

folder = "/home/sxn581/SRE_selection/"
setwd( folder )

# Construction of MR knots, and basis function -
source(paste( folder , "mr_knots.R", sep = "") )
source(paste( folder , "mr_radial_basis.R", sep = "") )
source(paste( folder , "mr_basis.R", sep = "") )

# TPS Basis -
source( paste( folder , "tps_basis.R", sep = "") )

# Likelihood Gradient and Fisher's Information for the SRE model -
sourceCpp(paste( folder , "Sigma_inv_det.cpp", sep = "") )
sourceCpp(paste( folder , "norm_cholesky.cpp", sep = "") )
sourceCpp(paste( folder , "like_grad_info_j.cpp", sep = "") )

# Armijo line search method -
sourceCpp(paste( folder , "BCGD_projection_j.cpp", sep = "") )
sourceCpp(paste( folder , "armijo_line_update.cpp", sep = "") )
source(paste( folder , "sre_grplasso_data.R", sep = "") )

# Scaling the coordinates to make a square domain :
# -------------------------------------------------

D0 = 2
n0 = nrow(coords)
s0 = ceiling(sqrt(n0))
Rep_overall = ncol(Y)
nres0 = 2
NC0 = 4
NCb0 = 0
Btype = "W"
Iter0 = 50

coords$x = s0*(coords$x - min(coords$x))/(max(coords$x)-min(coords$x))
coords$y = s0*(coords$y - min(coords$y))/(max(coords$y)-min(coords$y))

# Spatial knot locations, and Basis matrix :
# ------------------------------------------

# KInfo = mr_knots( D = D0 , coords , nres = nres0 , NC = NC0 )
# Rb = mr_basis( D = D0 , coords , knot_info = KInfo , basis_type = Btype, nres = nres0 )
# 
# OKnot_index = KInfo$Ordering # Ordered knot index from each resolutions
# nknots = KInfo$Kn # number of knots.

# LatticeKrig package -

ourDomain = apply( coords , 2 , "range" )
LKinfo = LKrigSetup( ourDomain , NC = NC0 , NC.buffer = NCb0 , nlevel = nres0 ,
                     alpha = c(1,0.5) , a.wght = 5.1 )
nknots = LKinfo$latticeInfo$m

Rb = as.matrix( LKrig.basis( x1 = coords, LKinfo = LKinfo , verbose = FALSE ) )

# TPS basis to compare with autoFRK package - 
# --------------------------------------------

# TPSinfo = tps_basis( D = 2 , coords , max_knots = n0 )
# nknots = TPSinfo$Kn
# 
# Rb = as.matrix( TPSinfo$Basis )

# Eigen values of the Basis matrix :
# ----------------------------------

Rbsq = crossprod.spam( Rb )
Eb_Rsq = eigen(Rbsq)$values
Eb_Rsq = Eb_Rsq[Eb_Rsq>0]
Cb1_n = min( Eb_Rsq )
Cb2_n = max( Eb_Rsq )
Kappab = Cb2_n/Cb1_n

# Centering the response :
# ------------------------

time = 1:Rep_overall
defr = 5 # 5 worked 10 did not work

# TB1 = bs( time , df = defr , Boundary.knots = c( 1 , ncol(Y) ) ) 
TB1 = bs( time , df = defr , intercept = TRUE , Boundary.knots = c( 1 , ncol(Y) ) ) 

# Exploratory DA:
# ---------------

# plot(time,TB1[,1])
# Add/Remove intercept here in the bs function
# Check alternative to the bs function 
# Dimension of TB1

# BETA = matrix( 0 , n0 , (defr+1) )
BETA = matrix( 0 , n0 , defr )
for( p in 1:n0 ){
  
  # p = 1
  
  # remove intercept from here
  
  # model = lm( Y[p,] ~ TB1 )
  model = lm( Y[p,] ~ TB1 - 1 )
  
  #plot(model,which=1)
  #plot( 1:40 , Y[p,] , type = "b")
  #lines( 1:40 , model$fitted.values , col="red" )
  
  BETA[p,] = model$coefficients
  
}

# smooth _ beta _bs over locations
# beta = BETA - t( matrix( apply( BETA , 2 , mean ) , (defr+1) , n0 ) )
beta = BETA - t( matrix( apply( BETA , 2 , mean ) , defr , n0 ) )

# E 
# E0 = Y - tcrossprod( beta , cbind(1,TB1) )
E0 = Y - tcrossprod( beta , TB1 )

Evec = as.vector( E0 )
tvec = rep( time , each = n0 )

model_sp = gam::gam( Evec ~ gam::s(tvec) )
mhat = matrix( model_sp$fitted.values , n0 , length(time) )

# smooth_beta = matrix( 0 , n0 , (defr+1) )
smooth_beta = matrix( 0 , n0 , defr )
for( q in 1:ncol(beta)){
  
  smooth_beta[,q] = gam::gam(beta[,q] ~ gam::s(coords$x , coords$y) )$fitted.values
  
}

# SP = Y - ( mhat + tcrossprod( smooth_beta , cbind(1,TB1) ) )
SP = Y - ( mhat + tcrossprod( smooth_beta , TB1 ) )

rm(model_sp)


for(i in 22:100){ sre_DOI_WR( N=i , E = SP , RB = Rb , RBCond = Kappab ) }



# Retrieve MSPE :
# ---------------

# NC = 2 , 3 , 4 , 5
# NC_buffer = 0 , 1

# Compare with the Lattice Krig prediction with our prediction for real data (and then with simulation)

D0 = 2
Rep0 = 75
# partition0 = c(5,6)

nknots = 65
# Btype = "LB"

M = 100
# Loss_F = rep( 0 , M )
# Loss_KL = rep( 0 , M )
rank_hat = rep( 0 , M )
mspe_AGL = rep( 0 , M )
mspe_LKRIG = rep( 0 , M )
rank_AFRK = rep( 0 , M )
mspe_AFRK = rep( 0 , M )
for( i in 1:M ){ 
  
  # i = 1
  
  # folder = "/home/sxn581/SRE_selection/MCMC_results_DOI/"
  
  fname = paste( folder , "DOI" , "Rep" , Rep0 , "Kn" , nknots , 
                 "MC" , i , '.RD' , sep = '' ) 
  
  load( fname ) 
  rank_hat[i] = output[1]
  mspe_AGL[i] = output[2]
  mspe_LKRIG[i] = output[3]
  rank_AFRK[i] = output[4]
  mspe_AFRK[i] = output[5]
  
}

( Avg_rank = mean( rank_hat ) )
( SD_rank = sd(rank_hat) )
( Avg_mspe_AGL = mean( mspe_AGL ) )
( SD_mspe_AGL = sd( mspe_AGL ) )
( Avg_mspe_LKRIG = mean( mspe_LKRIG ) )
( SD_mspe_LKRIG = sd( mspe_LKRIG ) )
( Avg_rank_AFRK = mean( rank_AFRK ) )
( SD_rank_AFRK = sd( rank_AFRK ) )
( Avg_mspe_AFRK = mean( mspe_AFRK ) )
( SD_mspe_AFRK = sd( mspe_AFRK ) )

#( Avg_LF = mean( Loss_F ) )
#( Avg_LKL = mean( Loss_KL ) )

