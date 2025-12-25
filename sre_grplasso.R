sre_grplasso = function( D , n , Rep , nres , NC , OKI , TKI , Kn , X , sigma , Phi_t , Rn , Xi , penalty_seq , Iter , stage ){
  
  # Arguments -
  # D is the dimension, where the locations are to be defined.
  # n is the number of observed spatial locations possibly irregular.
  # Rep is the number of independent copies of the spatial process.
  # nres is the number of resolution.
  # NC is the number of knots in coarsest resolution in one dimension.
  # OKI is the ordered knot index.
  # TKI is the true knot index.
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
  
  # D = D0 ; n = n0 ; Rep = Rep0 ; nres = nres0 ; NC = NC0 ; X = SP ; Kn = nknots ; Rn = Rt ; penalty_seq = tau1 ; Iter = Iter0
  # OKI = OKnot_index ; TKI = TKR_Index ; X = SP ; Kn = nknots ; sigma = sigmaT ; Phi_t = Phi_t_start ; Rn = Rt ; Xi = XI ; penalty_seq = tau1 ; stage = "G"
  
  O_index = 1:Kn
  Phi_init_t = Phi_t
  Phi0_t = Phi_t
  
  nCV = length( penalty_seq ) # number of cross validation.
  NLog_like = rep( 0 , nCV ) # collection of all cross validated Negative log likelihood.
  GIC = rep( 0 , nCV ) # collection of all cross validated GIC.
  rank_hat = rep( 0 , nCV ) # collection of all cross validated rank estimates.
  time = rep( 0 , nCV ) # collection of all computation times of cross validated iterations.
  nIteration = rep( 0 , nCV ) # collection of all cross validated number of iterations.
  TruePos = rep( 0 , nCV ) # collection of all cross validated true positives.
  FalsePos = rep( 0 , nCV ) # collection of all cross validated false positives.
  Cov_Est_acc2 = rep( 0 , nCV ) # collection of all cross validated overall covariance estimation Frobenius norm.
  # sigma_cv = rep( 0 , nCV ) # collection of all cross validated measurement error variance.
  Phi_cv = NULL # collection of all cross validated lower triangular parameter matrix.
  EIndex = NULL # collection of all estimated knot index.
  
  # k = 1
  
  for( k in 1:nCV ) {
    
    # k = k + 1
    
    time_k_start = Sys.time()
    
    # sigma_0 = sigma
    # sigma_0 = mean( apply( SP , 2 , var ) ) # initial measurement error variance.
    
    # Phi0_t = Phi_t
    # Phi_init_t = Phi_t
    
    i = 1
    acc.chck = 0.1
    while( ( acc.chck > 10^(-5) ) && ( i <= Iter ) ) {
      
      t1_i = Sys.time()
      
      Phim_t = Phi0_t # Counter to update the low rank Cholesky matrix.
      # sigma_m = sigma_0 # Counter to update the nugget variance.
      
      ProximalCGD = armijo_line_update( D , n , Kn , X , sigma , Phi_t = Phi0_t , 
                                        Phis_t = Phi_init_t , tau = penalty_seq[k] , 
                                        as.matrix( Rn ) , Rep , stage )
      
      Phi_init_t = Phim_t
      Phi0_t = ProximalCGD
      
      # If we estimate the nugget :
      # ---------------------------
      
      # Likelihood_info = neg_log_likelihood( n , X , sigma , Phi_t = Phi0_t , Kn , Rn , Rep )
      # Xi_0 = Likelihood_info$Cov
      # Xi_inv_0 = Likelihood_info$Cov_inv
      # SVC = Likelihood_info$SampleVC
      
      # D1 = crossprod.spam( Xi_inv_0 , Xi_0 - SVC )
      # gsigma = Rep * sum( D1 * Xi_inv_0 )
      # isigma = Rep * sum( Xi_inv_0 * Xi_inv_0 )
      
      # sigma_up = sigma_0 - gsigma/isigma
      # sigma_0 = sigma_up
      
      # ---------------------------
      
      acc.chck = sum( norm_cholesky( Phim_t - Phi0_t ) )
      
      i = i + 1
      
      t2_i = Sys.time() 
      ti = t2_i - t1_i
      
    } # end of iteration.
    
    Phi_hat = t(Phi0_t)
    lambda_Omega_hat = norm_cholesky( Phi0_t )
    EI = O_index[ which( lambda_Omega_hat > 0 ) ]
    
    Kn_t_res = function(t){(t*NC-(t-1))^D}
    
    CP = c( 0 , cumsum( Kn_t_res(1:nres) ) )
    EKR_Index = NULL
    TP = NULL
    FP = NULL
    for( q in 1:nres ){
      
      # q = 1
      
      EKR_Index[[q]] =  OKI[[q]][ EI[ which( CP[q] < EI & EI <= CP[q+1] ) ] - CP[q] ]
      TP[[q]] = length( intersect( TKI[[q]] , EKR_Index[[q]] ) )
      FP[[q]] = length( setdiff( EKR_Index[[q]] , intersect( TKI[[q]] , EKR_Index[[q]] ) ) )
      
    }
    
    EIndex[[k]] = EKR_Index
    
    TruePos[k] = do.call( 'sum' , TP )
    FalsePos[k] = do.call( 'sum' , FP )
    rank_hat[k] = TruePos[k] + FalsePos[k]
    
    A0 = tcrossprod.spam( Rn , Phi0_t ) # ( Rn*Phi )
    Xi_hat = sigma * diag.spam( n ) + tcrossprod.spam( A0 ) # ( sigma*In + Rn*Phi*t(Phi)*t(Rn) )

    NLog_like[k] = neg_log_likelihood( n , X , sigma , Phi_t = Phi0_t , Kn , as.matrix( Rn ) , Rep  )$NL
    GIC[k] = 2 * sum( EI ) + NLog_like[k]
    
    Phi_cv[[k]] = Phi_hat
    # sigma_cv[k] = sigma_0
    
    Cov_Est_acc2[k] = norm(Xi_hat - Xi,"F")
    
    time_k_end = Sys.time()
    time[k] = time_k_end - time_k_start
    nIteration[k] = i
    
  } # end of cross validation.
  
  Compare = data.frame( cbind( penalty_seq , time , nIteration , GIC , rank_hat , TruePos , FalsePos , Cov_Est_acc2 ) )

  min_GIC = which.min( GIC )
  TrueP = TruePos[ min_GIC ]
  FalseP = FalsePos[ min_GIC ]
  Phi_est = Phi_cv[[ min_GIC ]]
  ECV_Index = EIndex[[ min_GIC ]]
  # sigma_est = sigma_cv[ min_GIC ]
  A_hat = rank_hat[ min_GIC ]
  COV_EST_ACC2 = Cov_Est_acc2[ min_GIC ]
  
  return( list( Comp = Compare , 
                TrP = TrueP , 
                FaP = FalseP , 
                Phi_hat = Phi_est , 
                KR_hat = ECV_Index ,
                rn_hat = A_hat , 
                # COV_error1 = COV_EST_ACC1 , 
                COV_error2 = COV_EST_ACC2 ) )
  
}