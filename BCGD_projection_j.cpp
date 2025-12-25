#include <RcppArmadillo.h>
#include <Rcpp.h>

//#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

using namespace arma;

// [[Rcpp::export]]

List BCGD_projection_j( arma::mat Phi_t , arma::mat Phis_t , int Kn , double tau , arma::vec grad_vec , double info_max , int knot_id , char stage ){
  
  // Arguments -
  // Phi_t is the upper triangular matrix of the cholesky factor.
  // Phis_t is the starting value upper triangular matrix.
  // Kn is the number of knots.
  // tau is the global penalty parameter.
  // grad_vec is the gradient vector.
  // info_max is the maximum diagonal of the information matrix.
  // knot_id is the respective column of the cholesky matrix Phi.
  // stage is either "G" for first stage or "A" for adaptive stage.
    
  // Returns - 
  // the BC Gradient descent direction for a the $(knot_id)^th$ group or row of the lower triangular matrix Phi.
    
  // Values to test the function with -
  // Phi_t = Phi0_t ; Phis_t = Phi_init_t ; Kn = nknots ; tau = tau_seq[k] ; grad_vec = Gradient ; info_max = Fishers_Info_max ; knot_id = j
  
  arma::vec KKT_j = grad_vec - info_max * Phi_t.col(knot_id - 1) ;
  double KKTnorm_j = norm( KKT_j , 2 ) ;
  
  double psi_nj = 0 ;
  double tau_j = 0 ;
  
  if( stage == 'G' ) {

    psi_nj = sqrt( knot_id ) ;
    tau_j = tau * psi_nj ;
    // tau_j = tau ;

  } else {

    // psi_nj = Kn/( ( log(Kn) )*( Kn - knot_id + 1 ) ) ;
    // psi_nj = log( knot_id + 1 ) ; // WOR
    psi_nj = sqrt( knot_id ) ; // WR
    tau_j = tau * psi_nj / norm( Phis_t.col(knot_id - 1) , 2 ) ;
    
  }
  
  arma::vec dir_j( Kn , fill::zeros ) ;
  
  if( KKTnorm_j > tau_j ) { 
    
    dir_j = - ( 1 / info_max ) * ( grad_vec - ( tau_j * KKT_j / KKTnorm_j ) ) ;
    
  } else {
    
    dir_j = - Phi_t.col(knot_id - 1) ;
    
  }

  List result ;

  result["Dir"] = dir_j ;
  result["Tau"] = tau_j ;
  
  return result ;
  
}