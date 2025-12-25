true_knots = function( D , s , knot_info , nres , Model , rank ){
  
  # Arguments -
  # D is the dimension, where the locations are to  be defined.
  # s is the length of spatial grid in one direction. Default value is 1.
  # knot_info is the set of knot location and the ordering.
  # Model = 1 => true knots are subset of candidate knots.
  # Model = 2 => true knots are perturbation of a subset of candidate knots.
  # rank is the true rank of the low rank component.
  
  # Returns - 
  # The true knot locations (ordered).
  # The theta parameter for every resolution.
  # The diagonal distance.
  # The total number of knots.
  
  # Values to test the function with -
  # D = D0 ; s = sqrt(n) ; knot_info = KInfo ; nres = nres0 ; Model = TrueM ; rank = rn0
  
  KLoc = knot_info$Knot_Loc
  KIndex = knot_info$Ordering
  Ttheta = knot_info$Theta

  TKindex = NULL
  True_knot_loc = NULL
  True_theta = NULL
  for( t in 1:nres ){
    
    # t = 1
    KLt = KLoc[[t]]
    KIt = KIndex[[t]]

    if( rank%%2 == 0 ){ rnt = rank/2 }else{ rnt = (rank+1)/2 }
    if( nres == 1) {rnt = rank}
    
    set.seed( Model )
    # TKI = sample( KIt , rnt ) # If ordering is not used
    TKI = KIt[1:rnt] # If ordering is used
    
    index_match = which( match( as.numeric( row.names( KLt ) ) , TKI ) != "NA" )
    TKL = as.matrix( KLt[index_match,] )
    TKindex[[t]] = as.numeric( row.names( KLt )[ index_match ] )
    True_theta[t] = Ttheta[t]
    
    if( Model == 1 ){
      
      True_knot_loc[[t]] = TKL
      
    }
    if( Model == 2 ){
      
      # set.seed( rank + s )
      # Pb_x = runif( rank , -0.09 , 0.09 )
      # set.seed( 2*rank )
      # Pb_y = runif( rank , -0.09 , 0.09 )
      
      # Pb = cbind( Pb_x , Pb_y )
      # Tmp0 = do.call( 'rbind' , TKL )
      
      # True_knot_loc = Tmp0 + Pb
      
      # Tmp0 = TKL[[t]] # do.call( 'rbind' , TKL )
      Tmp1 = as.ppp( TKL , W = square(s) )
      Tmp2 = rjitter( Tmp1 , radius = 0.09 )
      
      True_knot_loc[[t]] = data.frame( cbind( Tmp2$x , Tmp2$y ) )
      
    }

  }
  
  return( list( Knot_Loc = True_knot_loc , 
                Theta = True_theta , 
                True_KRI = TKindex ) )
  
}