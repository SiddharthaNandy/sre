#include <RcppArmadillo.h>
#include <Rcpp.h>

//#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

using namespace arma;


// [[Rcpp::export]]


arma::vec norm_cholesky(arma::mat Phi_t){
  
  // Arguments -
  // Phi_t is an upper triangular matrix.
  
  // Returns - 
  // the 2-norm of every row or column of the positive semi-definite matrix corresponding to Phi_t.
  
  arma::uword Kn = Phi_t.n_cols ;
  arma::vec norms(Kn,fill::zeros) ;
  for (arma::uword i = 0; i < Kn; i++) {
    
    // norms(i) = norm( Phi_t.col(i) , 2 ) ;
    norms(i) = arma::norm(Phi_t.col(i),2);
    
  }
  
  return norms;
  
}