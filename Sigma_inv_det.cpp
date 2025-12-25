#include <RcppArmadillo.h>
#include <Rcpp.h>

//#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

using namespace arma;


// [[Rcpp::export]]
List Sigma_inv_det( int n , int Kn , double sigma , arma::mat Rn , arma::mat Phi_t ){
  
  arma::mat X = Rn*trans(Phi_t)/sqrt(sigma) ;
  arma::mat S0 = eye( n , n ) ;
  double logSD = 0.0 ;
  
  for( int i = 0 ; i < Kn ; ++i ){
    
    arma::vec Xi = X.col(i) ;
    arma::vec S1 = S0 * Xi ;
    
    double stemp = 1+dot( Xi , S1 ) ;
    
    S0 -= ( S1 * S1.t() )/stemp ;
    logSD += log(stemp) ;

  }
  
  S0 = S0/sigma ; 

  logSD += n*log(sigma) ;
    
  List result ;
  
  result["Inv"] = S0 ;
  result["LogDet"] = logSD ;
  
  return result ;

}