mr_basis = function( D , coords , knot_info , basis_type , nres ){
  
  # Arguments -
  # D is the dimension, where the locations are to  be defined.
  # coords is the location points.
  # knot_info contains the object mr_knots. 
  # basis_type is the type of basis matrix either Local_bisquare or Wendland.
  # nres is the number of resolutions using in MR.
  
  # Returnns -
  # the multi-resolution basis matrix for all resolutions using "mr_radial_basis.R". 

  # Values to test the function with -
  # D = D0 ; coords = sp_coords ; knot_info = KInfo ; basis_type = B_type ; nres = nres0
  
  RBM = NULL
  for( l in 1:nres ){
    
    # l = 1
    RB = mr_radial_basis( D , coords , knot_info , level = l , basis_type )$Rn
    RBM = cbind( RBM , RB )

  }
  
  return( Basis = as.spam( RBM , eps = 10^(-5) ) )

}