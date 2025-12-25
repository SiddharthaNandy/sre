mr_knots = function( D , coords , nres , NC ){
  
  # Arguments -
  # D is the dimension, where the locations are to  be defined.
  # coords is the observed spatial location points.
  # nres is the number of resolutions.
  # NC is the number of knots in coarsest resolution in one dimension.
  
  # Returns - 
  # The multi-resolution candidate knot locations (ordered).
  # The knot ordering.
  # The theta parameter for every resolution.
  # The diagonal distance.
  # The total number of knots.
  
  # Values to test the function with -
  # D = D0 ; coords = sp_coords ; nres = nres0 ; NC = NC0
  
  n = nrow( coords )
  s = sqrt( n )
  
  knot_loc = NULL # the knot locations or the centers.
  theta = NULL # the radius of the basis functions for each resolution.
  ordering = NULL
  OKL = NULL
  res_index = NULL
  Kn0 = 0
  
  for( t in 1:nres ) {

    # t = 1
    K0 = NULL
    for( d in 1:D ){
      
      # d = 1
      
      K0[[d]] = seq( 0 , s , length = t*NC - (t-1) )
      
    }
    
    knot_loc[[t]] = expand.grid( K0 )
    theta[t] = 1.5*min( dist( knot_loc[[t]] ) )
    
    # ordering[[t]] = 1:nrow(knot_loc[[t]])

    # distances = pointDistance( coords , knot_loc[[t]] , lonlat = F , allpairs = T )
    distances = crossdist( coords[,1] , coords[,2] , knot_loc[[t]][,1] , knot_loc[[t]][,2] )
    ordering[[t]] = order( apply( distances<theta[t] , 2 , sum ) , decreasing = T )
    
    OKL[[t]] = knot_loc[[t]][ordering[[t]], ]
    res_index[[t]] = rep( t , nrow( OKL[[t]] ) )
    
    Kn0 = Kn0 + nrow(knot_loc[[t]])
    
  }
  
  return( list( Knot_Loc = OKL , 
                Ordering = ordering , 
                Res_ind = res_index ,
                Theta = theta , 
                Kn = Kn0 ) )
  
}