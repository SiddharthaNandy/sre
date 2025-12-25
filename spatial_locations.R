spatial_locations = function( n , Irreg ){
  
  # Arguments - 
  # s is the length of spatial grid in one direction.
  # b is the boundary correction used to make sure the multi-resolution grid points contain all locations.
  # n is the number of observed spatial locations possibly irregular.
  # seed is used to permute the location.
  # Irreg is the amount of the perturbation.
  
  # Returns -
  # the simulated data locations
  
  # Values to test the function with -
  # n = n0 ; Irreg = Irreg0
  
  s = sqrt( n )
  
  set.seed( n )
  U = seq( 0 , s , length.out = s )
  
  # U = seq( b * s , ( 1 - b ) * s , length.out = sqrt(n) )
  
  S1 = expand.grid( U , U )
  
  if( Irreg == 0 ) {
    
    SL = S1
    
  }
  if( Irreg > 0 ) {
    
    # set.seed( n + Rep )
    # Purtb_x = runif( n , -Irreg , Irreg )
    # set.seed( 2*n + Rep )
    # Purtb_y = runif( n , -Irreg , Irreg )
    
    # Purtb = cbind( Purtb_x , Purtb_y )
    
    # SL = S1 + Purtb
    
    S2 = as.ppp( S1 , W = square(s) )
    S3 = rjitter( S2 , radius = Irreg )  
    
    SL = data.frame( cbind( S3$x , S3$y ) )

  }
  
  colnames(SL) = c("x","y")
  
  SL$x = sqrt(n)*(SL$x - min(SL$x))/(max(SL$x)-min(SL$x))
  SL$y = sqrt(n)*(SL$y - min(SL$y))/(max(SL$y)-min(SL$y))
  

  return(SL)
  
}