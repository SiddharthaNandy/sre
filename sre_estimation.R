sre_estimation = function( MC , D , n , Rep , coords , basis_type , nres , NC , NCb , Model , rank , Iter ) { 
  
  # Arguments - 
  # MC is the indexing number for simulations used to control the simulation seed number.
  # D is the dimension, where the locations are to be defined.
  # n is the number of observed spatial locations possibly irregular.
  # Rep is the number of independent copies of the spatial process.
  # coords is the observed spatial location points.
  # basis_type is the type of basis matrix either Wendland or local-bi-square.
  # nres is the number of resolution.
  # NC is the number of knots in coarsest resolution in one dimension.
  # Model indicates the number of true spatial random effects model of the process.
  # rank is the true rank of the low rank component.
  # Iter is the number of iterations.
  
  # Returns -
  # Cross validated Kriging estimate of the process, covariance matrix, and its rank.

  # Values to test the function with -
  # MC = N ; D = D0 ; n = n0 ;  Rep = Rep0 ; coords = sp_coords ; basis_type = B_type ; nres = nres0 ; NC = NC0 ; NCb = NCb0 ; Model = TrueM ; rank = rn0 ; Iter = Iter0
  
  # Knot info, and Basis matrix (Using knot ordering, but same scheme as Lattice-Kriging) :
  # ---------------------------------------------------------------------------------------
  
  KInfo = mr_knots( D , coords , nres , NC )
  Rt = mr_basis( D , coords , knot_info = KInfo , basis_type , nres )

  OKnot_index = KInfo$Ordering # Ordered knot index from each resolutions
  nknots = KInfo$Kn # number of knots.

  # Eigen values of the Basis matrix :
  # ----------------------------------
  
  Rsq = crossprod.spam( Rt )
  E_Rsq = eigen(Rsq)$values 
  E_Rsq = E_Rsq[E_Rsq>0]
  
  C1_n = min( E_Rsq )
  C2_n = max( E_Rsq )
  Kappa = C2_n/C1_n
  
  # C1_n ; C2_n ; Kappa
  
  # Data generation, basis matrix, and true parameter fixation :
  # ------------------------------------------------------------
  
  simulator = gen_data( MC , D , s = sqrt(n) , knot_info = KInfo , nres , coords , basis_type , Model , rank , n , Rep )
  
  SP = simulator$data # process.
  L_t = simulator$NLog_Like_Tr # likelihood value.
  PhiT_t = simulator$PhiTr_t # true upper triangular matrix of the cholesky factor.
  lambda_1 = max( norm_cholesky( as.matrix( PhiT_t ) ) ) # Largest column of PhiT_t
  sigmaT = simulator$sigmaTr_me # true measurement error.
  XI = simulator$true_cov # true overall covariance.
  XI_inv = simulator$true_cov_inv # inv of true overall covariance.
  TKR_Index = simulator$true_KR_index # true set
  RnT = simulator$Rntrue
  
  # Cholesky embedding :
  # --------------------

  hr = rank/2
  PhiT_t_Embd = cholesky_embedding( U = as.matrix(PhiT_t) , big_dim = nknots , TKR_Index,
                                    OKnot_index, j0 = hr )

  # Cross-validated estimation of the cholesky matrix in the First stage :
  # ----------------------------------------------------------------------

  set.seed( 1500 )
  # Lambda = (1:nknots)^(-1.5)
  Lambda = runif( nknots , 0.3 , 1.5 ) # (1:rank)^{-1.5}

  set.seed( 1500 * 2 )
  Omega0 = riwishart( nknots + 7 , diag( Lambda )  )
  Phi_t_start = chol(Omega0) # initial upper triangular parameter matrix.

  # Non-repetition -
  # C = 1 * ( 1/Kappa ) * seq( 100.0 , 1000.0 , length.out = 10 )
  # tau1 = C * sqrt ( sqrt( nknots ) )

  # Repetition -
  tau1 = (10^1) * ( 1/Kappa ) * seq( 100.0 , 1500.0 , length.out = 10 )
  # tau1 = C * sqrt( Rep * nknots )

  T1 = Sys.time()

  First_stage = sre_grplasso( D , n , Rep , nres , NC , OKI = OKnot_index , TKI = TKR_Index ,
                              Kn = nknots , X = SP , sigma = sigmaT , Phi_t = Phi_t_start ,
                              Rn = Rt , Xi = XI , penalty_seq = tau1 , Iter , stage = "G" )

  (T2 = Sys.time()-T1)

  Comp_GL = First_stage$Comp

  ## ERI_GL = First_stage$KR_hat

  TP_GL = First_stage$TrP
  FP_GL = First_stage$FaP

  Cov_Acc2_GL = First_stage$COV_error2/(C2_n*lambda_1*rank)

  # Cross-validated estimation of the cholesky matrix in the Adaptive stage :
  # -------------------------------------------------------------------------

  Phi_t_GL = t(First_stage$Phi_hat)

  # Non-repetition -
  # C = 10 * ( 1/Kappa ) * seq( 100.0 , 1000.0 , length.out = 10 )
  # tau2 = C * sqrt ( sqrt( nknots ) )

  # Repetition -
  C = (10^1) * ( 1/Kappa ) * seq( 10.0 , 100.0 , length.out = 10 )
  tau2 = C * sqrt( Rep * nknots )

  T1 = Sys.time()

  Adapt_stage = sre_grplasso( D , n , Rep , nres , NC , OKI = OKnot_index , TKI = TKR_Index ,
                              Kn = nknots , X = SP , sigma = sigmaT , Phi_t = Phi_t_GL ,
                              Rn = Rt , Xi = XI , penalty_seq = tau2 , Iter , stage = "A" )

  (T2 = Sys.time()-T1)

  Comp_AGL = Adapt_stage$Comp

  ## ERI_AGL = Adapt_stage$KR_hat

  TP_AGL = Adapt_stage$TrP
  FP_AGL = Adapt_stage$FaP

  Cov_Acc2_AGL = Adapt_stage$COV_error2/(C2_n*lambda_1*rank)

  res = PhiT_t_Embd - t(Adapt_stage$Phi_hat)

  return( list( TrueKnotIndex = TKR_Index ,
                BasisM_condN = Kappa ,
                # EI_GL = ERI_GL ,
                TrueP_GL = TP_GL ,
                FalseP_GL = FP_GL ,
                # FRnorm1_GL = Cov_Acc1_GL ,
                FRnorm2_GL = Cov_Acc2_GL ,
                # EI_AGL = ERI_AGL ,
                TrueP_AGL = TP_AGL ,
                FalseP_AGL = FP_AGL ,
                # FRnorm1_AGL = Cov_Acc1_AGL ,
                FRnorm2_AGL = Cov_Acc2_AGL ,
                res = res ) ) # ,
                # z_score = z ) ) #,
  
}