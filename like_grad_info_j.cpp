#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Sigma_inv_det.cpp>

//#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

using namespace arma;


// [[Rcpp::export]]
List grad_info_j( int n , arma::mat X , double sigma , arma::mat Phi_t , int Kn , arma::mat Rn , int Rep , int knot_id ){
  
  // Arguments -
  // n is the number of observed spatial locations possibly irregular.
  // X is the observed spatial process.
  // sigma is the nugget value.
  // Phi_t is the upper triangular matrix as a spam object.
  // Kn is the number of knots.
  // Rn is the basis matrix as a spam object.
  // Rep is the number of independent repetition. 
  // knot_id is the respective column of the cholesky matrix Phi.
  // We use the structure, Delta_Omega_ij = A + A', where A is a matrix with only one non-zero column.
  
  // Returns -
  // The gradient vector of the likelihood with respect to the (knot_id)^th column.
  // The largest diagonal of the information matrix.
  
  // Values to test the function with -
  // n = n0 ; X = SP ; sigma = sigmaT ; Phi_t = Phi_start_t ; Kn = nknots ; Rn = Rt ; Rep = Rep0 ; knot_id = j
  
  List Sig_ID = Sigma_inv_det( n , Kn , sigma , Rn , Phi_t ) ;
  arma::mat Xi_inv = Sig_ID["Inv"] ;
  double logdet = Sig_ID["LogDet"] ;
  
  arma::mat Xi_inv_sqrt = sqrtmat_sympd( Xi_inv ) ;
  arma::mat SVC = X*trans(X)/Rep ;
  double NLogLike = Rep * accu( Xi_inv%SVC ) + Rep * logdet ;
  
  arma::mat T3 = SVC*Xi_inv_sqrt ; 
  arma::mat T4 = eye(n,n) - Xi_inv_sqrt*T3 ; 
  arma::mat T5 = Xi_inv_sqrt*Rn ;
  
  // arma::vec T5_j_col = T5*Phi_t.col(knot_id-1) ; 
  
  arma::rowvec T5jr = trans( T5.col(knot_id-1) ) ;
  arma::mat Phi = trans(Phi_t) ;
  arma::vec T5a(n,fill::zeros) ;
  
  arma::vec T5_i_row( Kn , fill::zeros ) ;
  arma::vec gtemp( Kn , fill::zeros) ;
  arma::vec itemp( Kn , fill::zeros) ;
  arma::mat T6(n,n,fill::zeros) ;
  arma::mat D_ij(n,n,fill::zeros) ;
  
  for(int i=0; i<knot_id;++i){
    
    T5a = T5*Phi.col(i) ;
    T6 = T5a*T5jr;
    
    // T5_i_row = T5.col(i) ;
    // T6 = T5_j_col*T5_i_row.t() ;
    gtemp(i) = Rep * 2 * accu( T4%T6 ) ;
    
    D_ij = T6 + trans(T6);
    itemp(i) = Rep * accu( D_ij%D_ij ) ;
    
  } // # for.
  
  arma::vec I0(1,fill::zeros) ;
  arma::vec I1(1,fill::zeros) ;
  I0 = pow(10,-4) ;
  I1 = pow(10,8) ;
  
  arma::vec I2(1,fill::zeros)  ;
  I2 = max( join_cols( itemp , I0 ) ) ;
  
  arma::vec I3 = join_cols( I2 , I1 ) ;
  
  double info_max = min( I3 ) ;
  
  List result;
  
  result["NL"] = NLogLike ;
  result["grad"] = gtemp ;
  result["info"] = info_max ;
  result["FI"] = itemp ;
  
  return result ;
  
}



// [[Rcpp::export]]
List neg_log_likelihood( int n , arma::mat X , double sigma , arma::mat Phi_t , int Kn , arma::mat Rn , int Rep ){
  
  // Arguments -
  // n is the number of observed spatial locations possibly irregular.
  // X is the observed spatial process.
  // sigma is the nugget value.
  // Phi_t is the upper triangular matrix as a spam object.
  // Kn is the number of knots.
  // Rn is the basis matrix as a spam object.
  // Rep is the number of independent repetition. 
  
  // Returns -
  // The gradient vector of the likelihood with respect to the (knot_id)^th column.
  // The largest diagonal of the information matrix.
  
  // Values to test the function with -
  // n = n0 ; X = SP ; sigma = sigmaT ; Phi_t = Phi_start_t ; Kn = nknots ; Rn = Rt ; Rep = Rep0 ; knot_id = j
  
  List Sig_ID = Sigma_inv_det( n , Kn , sigma , Rn , Phi_t ) ;
  arma::mat Xi_inv = Sig_ID["Inv"] ;
  double logdet = Sig_ID["LogDet"] ;
  
  arma::mat SVC = X*trans(X)/Rep ;
  double NLogLike = Rep * accu( Xi_inv%SVC ) + Rep * logdet ;
  
  List result ;
  
  result["NL"] = NLogLike ;
  result["SampleVC"] = SVC ;
  // result["Cov"] = Xi ;
  result["Cov_inv"] = Xi_inv ;
  
  return result ;
  
  // return NLogLike ;
  
}