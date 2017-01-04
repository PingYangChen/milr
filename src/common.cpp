#include "common.h"

// function to check whether the input data with correct type
void chk_mat(const mat& x, const std::string& varName){
  if (!is_finite(x))
    Rcpp::stop(varName + " must be numerical.\n");
}

//' logit link function
//'
//' calculate the values of logit link
//' 
//' @param X A matrix, the design matrix.
//' @param beta A vector, the coefficients.
//' @return An vector of the values of logit link.
// [[Rcpp::export]]
arma::vec logit(const arma::mat& X, const arma::vec& beta){
  chk_mat(X, "X");
  chk_mat(beta, "beta");
  return pow(1.0 + exp(-X * beta), -1.0);
}
