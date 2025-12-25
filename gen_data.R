gen_data = function( MC , D , s , knot_info , nres , Crds , basis_type , Model , rank , n , Rep ){
  
  # Arguments - 
  # MC is the indexing number for simulations used to control the simulation seed number.
  # D is the dimension, where the locations are to be defined.
  # s is the length of spatial grid in one direction.
  # Crds is the location points.
  # Model value is the true spatial random effects model of the process.
  # n is the number of observed spatial locations possibly irregular.
  # Rep is the number of independent repetition. 
  
  # Returns -
  # Generated data for nd-MC simulation.
  
  # Values to test the function with -
  # MC = N ; D = D0 ; s = sqrt(n0) ; knot_info = KInfo ; nres = nres0 ; Crds = coords ; basis_type = B_type ; Model = TrueM ; rank = rn0 ; n = n0 ; Rep = Rep0
  # knot_info = KInfo ; Crds = coords
  
  True_Kn = true_knots( D , s , knot_info , nres , Model , rank )
  true_knot_index = True_Kn$True_KRI
  
  set.seed( 1500 )
  # lambda = (1:rank)^(-2.5)
  # lambda = sort( runif( rank , 7 , 14 ) , decreasing = T ) 
  # lambda = sort( runif( rank , 1.5 , 2.5 ) , decreasing = T ) # This is our choice used for repetition case
  lambda = sort( runif( rank , 0.5 , 1.5 ) , decreasing = T ) # This choice is more supportive of adaptive selection
  
  set.seed( 1500 * 2 )
  Omega = riwishart( rank + 10 , diag( lambda )  )
  
  PhiTrue_t = chol(Omega)
  
  PhiTrue_t = as.spam( PhiTrue_t )
  Rnt = mr_basis( D , Crds , knot_info = True_Kn , basis_type , nres )
  # True_basis = mr_radial_basis( D , Crds , knot_info = True_Kn , level = "N" , basis_type )
  # Rnt = True_basis
  
  M0 = tcrossprod.spam( Rnt , PhiTrue_t ) # ( Rn*Phi )
  M1 = crossprod.spam( M0 ) # ( t(Phi)*t(Rn)*Rn*Phi )
  
  # choice for SNR: 
  sigma_me = sum( diag( tcrossprod.spam( M0 ) ) )/n # or, user defined choices : 15 , 3 , 0.5 , 0.1 , 0.05
  
  M2 = sigma_me * diag.spam( rank ) + M1 # ( sigma*I_Kn + t(Phi)*t(Rn)*Rn*Phi )
  M3 = t(M0)
  M4 = crossprod.spam( M3 , solve.spam( M2 , M3 ) ) # ( Rn*Phi*( sigma*I_Kn + t(Phi)*t(Rn)*Rn*Phi )^(-1)*t(Rn*Phi) )
  
  Xi = sigma_me * diag.spam(n) + tcrossprod.spam( M0 ) # ( sigma*In + Rn*Phi*t(Phi)*t(Rn) )
  Xi_inv = ( 1/sigma_me ) * ( diag.spam(n) - M4  ) # (1/sigma)*( I_n - Rn*Phi*( sigma*I_Kn + t(Phi)*t(Rn)*Rn*Phi )^(-1)*t(Phi)*t(Rn) )
  
  log_det = determinant.spam(Xi)$modulus # log( det( sigma*In + Rn*Phi*t(Phi)*t(Rn) ) )
  
  set.seed(MC)
  P = t( as.matrix( rmvnorm( Rep , rep( 0 , n ) , diag( sigma_me , n ) ) ) )
  set.seed(2*MC+3)
  E = t( as.matrix( rmvnorm( Rep , rep( 0 , rank ) , diag( rank ) ) ) )
  E = crossprod.spam( t(M0) , E )
  Y = as.matrix( E + P )
  
  SVarCov = tcrossprod.spam( Y ) / Rep
  Nloglikelihood = as.vector( Rep * tr( crossprod.spam( Xi_inv , SVarCov ) ) + Rep * log_det )
  
  return( list( data = Y ,
                NLog_Like_Tr = Nloglikelihood , 
                true_cov = Xi , 
                true_cov_inv = Xi_inv ,
                PhiTr_t = PhiTrue_t , 
                sigmaTr_me = sigma_me , 
                true_KR_index = true_knot_index , 
                Rntrue = Rnt ) )
  
}