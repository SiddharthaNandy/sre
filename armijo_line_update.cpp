#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <like_grad_info_j.cpp>
#include <BCGD_projection_j.cpp>
#include <norm_cholesky.cpp>

// #ifdef _OPENMP
// #include <omp.h>
// #endif

//// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

using namespace arma;


// [[Rcpp::export]]
arma::mat armijo_line_update( int D , int n , int Kn , arma::mat X , double sigma , arma::mat Phi_t , arma::mat Phis_t , double tau , arma::mat Rn , int Rep , char stage ){
    
  // Arguments -
  // D is the dimension, where the locations are to be defined.
  // Kn is the total number of knots.
  // X is the observed spatial process.
  // sigma is the nugget parameter.
  // Phi_t is an upper triangular matrix.
  // Phis_t is the starting value upper triangular matrix.
  // tau is the global penalty parameter.
  // Rn is the basis matrix as a spam object.
  // Rep is the number of independent copies of the spatial process.
  // stage is either "G" for first stage or "A" for adaptive stage.

  // Returns -
  // armijo line search update.
    
  // Values to test the function with -
  // D = D0 ; n = n0 ; knot_info = KInfo ; X = SP ; sigma = sigmaT ; Phi_t = Phi_start_t ; tau = tau1(1) ; Rn = Rt ; Rep = Rep0 ; stage = "G"
  // Phi_t = Phi0_t ; Phis_t = Phi_init_t ; tau = penalty_seq(1) ; Rn = as.matrix( Rt ) ; Rep = Rep0 ; stage = "G"

  int knot_id = 0 ;
  
  // arma::mat PhiDirection_t( Kn , Kn , fill::zeros ) ;
  // arma::mat PhiD_t( Kn , Kn , fill::zeros ) ;
  
  double gamma = 0 ; // (0<=gamma<1)
  double Delta_j = 0 ;
  
  double alpha_j = 1 ;
  
  double Qdiff_j = 0 ;
  int counter = 0 ;
  
  double NLL1 = 0 ;
  double Pen1 = 0 ;
  double NLL2 = 0 ;
  double Pen2 = 0 ;
  
  double Q1 = 0 ;
  double Q2 = 0 ;
  
  // #pragma omp parallel for if(_OPENMP)
  for(int j=0 ; j<Kn ; ++j ){
    
    knot_id = j + 1 ;
    
    List GI_j = grad_info_j( n , X , sigma , Phi_t , Kn , Rn , Rep , knot_id ) ;
    arma::vec grad_vec = GI_j["grad"] ;
    double info_max = GI_j["info"] ;

    List BCGD_Proj = BCGD_projection_j( Phi_t , Phis_t , Kn , tau , grad_vec , info_max , knot_id , stage ) ;
    arma::vec BCGD_direction = BCGD_Proj["Dir"] ;
    double PENALTY_j = BCGD_Proj["Tau"] ;

    if( norm( BCGD_direction , 2 ) == 0 ){ 
      
      // NOT updating the upper triangular matrix using Armijo line search:
      Phi_t.col(j) = BCGD_direction ; 
      
    } else {
      
      // Updating the upper triangular matrix using Armijo line search:
      arma::mat PhiDirection_t( Kn , Kn , fill::zeros ) ;
      // arma::mat PhiD_t( Kn , Kn , fill::zeros ) ;
      
      PhiDirection_t.col(j) = BCGD_direction ;
      
      Delta_j = dot( BCGD_direction , grad_vec ) + gamma * info_max * dot( grad_vec , grad_vec ) 
              + PENALTY_j * (  norm( Phi_t.col(j) + BCGD_direction , 2 ) - norm( Phi_t.col(j) , 2 ) ) ;
      
      // For Simulation study I used : step size = 10^(-(NC^nres))
      // For Real data I used : step size = 10^(-(NC^(nres+2))) 

      while ( Qdiff_j > alpha_j*(std::pow(10,-9))*Delta_j ) {
        
        counter = counter + 1 ;
        
        arma::mat PhiD_t = Phi_t + alpha_j * PhiDirection_t ;
        
        List NL_one = neg_log_likelihood( n , X , sigma , Phi_t , Kn , Rn , Rep ) ;
        
        NLL1 = NL_one["NL"] ; // NL_one["NL"] ; // GI_j["NL"] ;   
        Pen1 = PENALTY_j * sum( norm_cholesky( Phi_t ) ) ;
        
        List NL_two = neg_log_likelihood( n , X , sigma , PhiD_t , Kn , Rn , Rep ) ;
        
        NLL2 = NL_two["NL"] ;
        Pen2 = PENALTY_j * sum( norm_cholesky( PhiD_t ) ) ;
        
        Q1 = NLL1 + Pen1 ;
        Q2 = NLL2 + Pen2 ;
        
        Qdiff_j = Q2 - Q1 ;
        alpha_j = 0.5 * alpha_j ;
        
      } // while
      
      // updating the upper triangular matrix using Armijo line search:
      Phi_t.col(j) = Phi_t.col(j) + 2 * alpha_j * BCGD_direction ; 

    } // if_else.
    
  } // for.
  
  return Phi_t ;
  
} 