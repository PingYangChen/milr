#include "common.h"

// [[Rcpp::export]]
double softmaxlogL(const arma::vec& bag, const arma::mat& X, const arma::vec& Z, 
                   const arma::vec& beta, const double& alpha){
  arma::vec p_vec = logit(X, beta);
  p_vec.elem(find(p_vec < 1e-7)).fill(1e-7);
  p_vec.elem(find(p_vec > 1-1e-7)).fill(1-1e-7);
  arma::vec p_bag = arma::zeros<arma::vec>(Z.n_elem);
  arma::vec denominator = arma::zeros<arma::vec>(Z.n_elem);
  double tmp;
  
  for (arma::uword i = 0; i < bag.n_elem; i++) {
    tmp = std::exp(alpha * p_vec[i]);
    p_bag[bag[i] - 1] += p_vec[i] * tmp;
    denominator[bag[i] - 1] += tmp;
  }
  p_bag /= denominator;
  
  double sum1 = -arma::sum(Z % arma::log(p_bag) + (1 - Z) % arma::log(1 - p_bag));
  return sum1;
}
