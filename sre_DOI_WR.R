sre_DOI_WR = function( N , E , RB , RBCond ) {
  
  # rm(list=ls())
  
  # N = 91 ; E = SP ; RB = Rb ; RBCond = Kappab # Comment this out to run in parallel
  
  T1 = Sys.time()
  
  set.seed(N)
  Rep_train = 75
  Rep_test = Rep_overall - Rep_train
  
  train_set = sort( sample.int( Rep_overall , Rep_train ) )
  E_train = E[ , train_set ]
  Scov_train = tcrossprod( E_train )/Rep_train
  
  test_set = setdiff( 1:Rep_overall , train_set )
  E_test = E[ , test_set ]
  Scov_test = tcrossprod( E_test )/Rep_test
  
  # Setting the starting value of Cholesky and Covariance matrix:
  # -------------------------------------------------------------
  
  # TV = mean( apply( E_train , 2 , var) ) # var(E) # Total variability of the process.
  TV = sum(diag(Scov_train))/n0 # var(E) # Total variability of the process.
  sf = 1 # Increasing the factor as Kn increases shows better and faster accuracy.
  
  sigmaT = sf*TV # nugget starting value
  
  set.seed( 1500 )
  Lambda = runif( nknots , 0.3 , 1.5 ) # (1:rank)^{-1.5}
  
  set.seed( 1500 * 2 )
  Omega0 = riwishart( nknots + 9 , diag( Lambda )  )
  Phi_start_t = chol(Omega0) # initial upper triangular parameter matrix.
  
  # Cross-validated estimation of the cholesky matrix in the First stage :
  # ----------------------------------------------------------------------
  
  # MR knots LatticeKrig
  # NC0 = 4 ; (10^(2*NC0-1)) ; (10^(2*NC0-2)) 
  # NC0 = 5 ; (10^(2*NC0+1)) ; (10^(2*NC0))
  # NC0 = 6 ; (10^(2*NC0+5)) ; (10^(2*NC0+4))
  # NC0 = 7 ; (10^(2*NC0+6)) ; (10^(2*NC0+5))
  
  # check the order of the tuning  parameter
  
  tau1 = (10^(2*NC0-1)) * ( 1/RBCond ) * seq( 500.0 , 5000.0 , length.out = 10 ) 
  # tau1 = C * sqrt( Rep_train * nknots )
  
  # TPS
  # tau1 = (10^15) * ( 1/RBCond ) * seq( 50.0 , 2000.0 , length.out = 10 ) 
  
  T1 = Sys.time()
  
  First_stage = sre_grplasso_data( D = D0 , n = n0 , Rep = Rep_train , nres = nres0 , 
                                   Kn = nknots , X = E_train , sigma = sigmaT , 
                                   Phi_t = Phi_start_t , Rn = as.matrix(RB) , penalty_seq = tau1 , 
                                   Iter = Iter0 , stage = "G" )
  
  T2 = Sys.time()
  (Computation.time = T2 - T1)
  
  Comp_GL = First_stage$Comp
  ERI_GL = First_stage$KR_hat
  
  Phi_f_t = t(First_stage$Phi_hat)
  sigma_f = First_stage$sigma_hat
  
  Btemp = tcrossprod.spam( RB , Phi_f_t ) # ( Rn*Phi )
  Xi_GL = sigma_f * diag.spam( n0 ) + tcrossprod.spam( Btemp ) # ( sigma*In + Rn*Phi*t(Phi)*t(Rn) )
  
  # Cross-validated estimation of the cholesky matrix in the Adaptive stage :
  # -------------------------------------------------------------------------
  
  Phi_adapt_t = t(First_stage$Phi_hat)
  
  # MR Latticekriging - 
  C = (10^(2*NC0-2)) * ( 1/RBCond ) * seq( 500.0 , 5000.0 , length.out = 10 ) 
  tau2 = C * sqrt( Rep_train * nknots )
  
  # TPS
  # C = (10^13) * ( 1/RBCond ) * seq( 50.0 , 2000.0 , length.out = 10 )
  # tau2 = C * sqrt( Rep_train * nknots )
  
  T3 = Sys.time()
  
  Adapt_stage = sre_grplasso_data( D = D0 , n = n0 , Rep = Rep_train , nres = nres0 , 
                                   Kn = nknots , X = E_train , sigma = sigma_f , 
                                   Phi_t = Phi_f_t , Rn = as.matrix(RB) , penalty_seq = tau2 , 
                                   Iter = Iter0 , stage = "A" )
  
  T4 = Sys.time()
  (Computation.time = T4 - T3)
  
  Comp_AGL = Adapt_stage$Comp
  ERI_AGL = Adapt_stage$KR_hat
  rank_AGL = Adapt_stage$rn_hat
  
  Phi_s_t = t(Adapt_stage$Phi_hat)
  sigma_s = Adapt_stage$sigma_hat
  
  Xi_AGL_inv = Sigma_inv_det( n = n0 , Kn = nknots , sigma = sigma_s , Rn = as.matrix(RB) , 
                          Phi_t = Phi_s_t )$Inv
  
  # Predicting using 'LatticeKrig'
  # ------------------------------
  
  LKobj = LKrig( coords , E_test , LKinfo = LKinfo , lambda = 0.1 )
  
  # Predicting using 'AutoFRK'
  # --------------------------
  
  AFRK = autoFRK( Data = E_train , loc = coords , mu = 0, D = diag.spam(NROW(E_train)),
                  G = NULL, finescale = FALSE, maxit = 50, tolerance = 0.1^6, maxK = NULL ,
                  Kseq = NULL, method = "EM", maxknot = 5000 )

  FMFt = tcrossprod.spam( tcrossprod.spam( AFRK$G , AFRK$M ) , AFRK$G )
  Xi_AFRK = AFRK$s * diag(n0) + FMFt
  Xi_AFRK_inv = solve( Xi_AFRK )
  rank_AFRK = rankMatrix( AFRK$M )
  
  # MSPE :
  # ------
  
  Ehat_test_AGL = E_test - Xi_AGL_inv %*% E_test
  Ehat_test_LKRIG = predict( LKobj , xnew = LKobj$x )
  Ehat_test_AFRK = tcrossprod( FMFt , Xi_AFRK_inv ) %*% E_test
  
  MSPE_AGL = sum((Ehat_test_AGL - E_test)^2)/(n0*Rep_test)
  MSPE_LKRIG = sum((Ehat_test_LKRIG - E_test)^2)/(n0*Rep_test)
  MSPE_AFRK = sum((Ehat_test_AFRK - E_test)^2)/(n0*Rep_test)
  
  # output = c( rank_AGL , MSPE_AGL )
  # output = c( rank_AGL , MSPE_AGL , MSPE_LKRIG )
  output = c( rank_AGL , MSPE_AGL , MSPE_LKRIG , rank_AFRK , MSPE_AFRK )
  
  fname = paste( folder , "DOI" , "Rep" , Rep_train , "Kn" , nknots , "MC" , N , '.RD' , sep = '' ) 
  save( output , file = fname ) 
  
  T2 = Sys.time()
  Computation.time = T2 - T1
  
  print(N)
  print(Computation.time)
  
  # return( list ( LF = LossF , LKL = LossKL ) )

}