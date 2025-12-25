mr_radial_basis = function( D , coords , knot_info , level , basis_type ){
  
  # Arguments -
  # D is the dimension, where the locations are to  be defined.
  # coords are the set generate location points.
  # knot_info contains the object mr_knots. 
  # level is the resolution level.
  # basis_type is the type of the basis funciton used for the estimation.
  
  # Returns -
  # basis matrix of a particular level of resolution.
  
  # Values to test the function with -
  # D = D0 ; coords = coords ; knot_info = KInfo ; level = "N" ; basis_type = B_type
  
  n = nrow( coords )
  
  if( level == "N" ){ 
    
    Knots = knot_info$TKL
    radius = knot_info$TT
    ordering = NULL
    
  }
  else{ 

    Knots = knot_info$Knot_Loc[[level]]
    radius = knot_info$Theta[level]
    ordering = knot_info$Ordering[[level]]
    
  }
    
  # B1 = pointDistance( coords , Knots , lonlat = F , allpairs = T ) / radius
  B1 = crossdist( coords[,1] , coords[,2] , Knots[,1] , Knots[,2] ) / radius
  B1[ B1 > 1 ] = 0
  
  if ( basis_type == "LB" ) {
    
    B2 = ( 1 - ( B1 )^2 )^2
    
  }
  if ( basis_type == "W" ) {
    
    B2 = ( ( 1 - B1 )^6 )*( 35*B1^2 + 18*B1 + 3 ) / 3
    
  }
  
  B2[ B2 == 1 ] = 0
  mr_RB = B2

  # return( Rn = mr_RB ) 
  return( list( Rn = mr_RB , OK_Index = ordering ) )
  
}