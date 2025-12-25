sre_grplasso_data = function( D , n , Rep , nres , Kn , X , sigma , Phi_t , Rn , penalty_seq , Iter , stage ){
  
  # Arguments -
  # D is the dimension, where the locations are to be defined.
  # n is the number of observed spatial locations possibly irregular.
  # Rep is the number of independent copies of the spatial process.
  # nres is the number of resolution.
  # knot_info contains knot information either from radial_basis ( for no resolution ) or mr_knots ( for multiple resolutions ). 
  # Kn is the total number of knots.
  # X is the observed spatial process.
  # sigma is the nugget parameter.
  # Phi_t is an upper triangular matrix as a spam object.
  # Rn is the basis matrix as a spam object.
  # penalty_seq is the sequence of penalty parameter.
  # Iter is the number of iteration.
  # stage is either "G" for first stage or "A" for adaptive stage.
  
  # Returns -
  # group LASSO estimate.
  
  # D = D0 ; n = n0 ; Rep = Rep_train ; nres = nres0 ;  
  # Kn = nknots ; X = E_train ; sigma = sigmaT ; Phi_t = Phi_start_t ; 
  # Rn = RB ; penalty_seq = tau1 ; Iter = Iter0 ; stage = "G"
  
  O_index = 1:Kn
  # Phi_init_t = Phi_t
  
  nCV = length( penalty_seq ) # number of cross validation.
  NLog_like = rep( 0 , nCV ) # collection of all cross validated Negative log likelihood.
  GIC = rep( 0 , nCV ) # collection of all cross validated GIC.
  rank_hat = rep( 0 , nCV ) # collection of all cross validated rank estimates.
  time = rep( 0 , nCV ) # collection of all computation times of cross validated iterations.
  nIteration = rep( 0 , nCV ) # collection of all cross validated number of iterations.
  sigma_cv = rep( 0 , nCV ) # collection of all cross validated measurement error variance.
  Phi_cv = NULL # collection of all cross validated lower triangular parameter matrix.
  EIndex = NULL # collection of all estimated knot index.
  
  for( k in 1:nCV ) {
    
    # k = 1
    
    time_k_start = Sys.time()
    
    sigma_0 = sigma
    # sigma_0 = mean( apply( SP , 2 , var ) ) # initial measurement error variance.
    
    Phi0_t = Phi_t
    Phi_init_t = Phi_t
    
    i = 1
    acc.chck = 0.1
    while( ( acc.chck > 10^(-5) ) && ( i <= Iter ) ) {
      
      t1_i = Sys.time()
      
      Phim_t = Phi0_t # Counter to update the low rank cholesky matrix.
      sigma_m = sigma_0 # Counter to update the nugget variance.
      
      ProximalCGD = armijo_line_update( D , n , Kn , as.matrix(X) , sigma = sigma_0 , Phi_t = Phi0_t , 
                                        Phis_t = Phi_init_t , tau = penalty_seq[k] , Rn , Rep , stage )
      
      Phi_init_t = Phim_t
      Phi0_t = ProximalCGD
      
      # If we estimate the nugget :
      # ---------------------------
      
      Likelihood_info = neg_log_likelihood( n , as.matrix(X) , sigma = sigma_0 , Phi_t = Phi0_t , Kn , Rn , Rep )
      # Xi_0 = Likelihood_info$Cov
      Xi_inv_0 = Likelihood_info$Cov_inv
      SVC = Likelihood_info$SampleVC
      
      D1 = diag(n) - crossprod( Xi_inv_0 , SVC )
      gsigma = Rep * sum( D1 * Xi_inv_0 )
      isigma = Rep * sum( Xi_inv_0 * Xi_inv_0 )
      
      sigma_up = sigma_0 - gsigma/isigma
      sigma_0 = sigma_up
      
      # ---------------------------
      
      acc.chck = abs( sigma_m - sigma_0 ) + sum( sqrt( norm_cholesky( Phim_t - Phi0_t ) ) )
      
      i = i + 1
      
      t2_i = Sys.time() 
      ti = t2_i - t1_i
      # print( c( ti , acc.chck ) )
      
    } # end of iteration.
    
    Phi_hat = t(Phi0_t)
    lambda_Omega_hat = norm_cholesky( Phi0_t )
    # lambda_Omega_hat = round( norm_cholesky( Phi0_t ) )
    EI = O_index[ which( lambda_Omega_hat > 0 ) ]
    
    # alternative way to EKR_Index ( we do not need TP and FP)
    
    rank_hat[k] = length( EI )

    A0 = tcrossprod.spam( Rn , Phi0_t ) # ( Rn*Phi )

    Omega_hat = tcrossprod.spam( Phi_hat )
    Xi_hat = sigma_0 * diag.spam( n ) + tcrossprod.spam( A0 ) # ( sigma*In + Rn*Phi*t(Phi)*t(Rn) )
    # Xi_inv_hat = ( 1/sigma_0 ) * ( diag.spam( n ) - A3  )
    
    NLog_like[k] = neg_log_likelihood( n , X , sigma = sigma_0 , Phi_t = Phi0_t , Kn , Rn , Rep  )$NL
    GIC[k] = 2 * sum( EI ) + NLog_like[k]

    Phi_cv[[k]] = Phi_hat
    sigma_cv[k] = sigma_0
    

    time_k_end = Sys.time()
    time[k] = time_k_end - time_k_start
    nIteration[k] = i
    
  } # end of cross validation.
  
  Compare = data.frame( cbind( penalty_seq , time , nIteration , GIC , rank_hat ) )

  min_GIC = which.min( GIC )
  Phi_est = Phi_cv[[ min_GIC ]]
  ECV_Index = EIndex[[ min_GIC ]]
  sigma_est = sigma_cv[ min_GIC ]
  A_hat = rank_hat[ min_GIC ]

  return( list( Comp = Compare , 
                # TrP = TrueP , 
                # FaP = FalseP , 
                sigma_hat = sigma_est ,
                Phi_hat = Phi_est , 
                KR_hat = ECV_Index ,
                rn_hat = A_hat ) ) #, 
  # COV_error1 = COV_EST_ACC1 , 
  # COV_error2 = COV_EST_ACC2 ) )
  
}