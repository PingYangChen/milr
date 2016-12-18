#include "common.h"

// [[Rcpp::export]]
double softmaxlogL(const arma::vec& bag, const arma::mat& X, const arma::vec& Z, 
                   const arma::vec& beta, const double& alpha){
  vec p_vec = logit(X, beta);
  p_vec.elem(find(p_vec < 1e-7)).fill(1e-7);
  p_vec.elem(find(p_vec > 1-1e-7)).fill(1-1e-7);
  vec p_bag = zeros<vec>(Z.n_elem);
  vec denominator = zeros<vec>(Z.n_elem);
  double tmp;
  
  for(uword i = 0; i < bag.n_elem; i++)
  {
    tmp = exp(alpha * p_vec[i]);
    p_bag[bag[i] - 1] += p_vec[i] * tmp;
    denominator[bag[i] - 1] += tmp;
  }
  p_bag /= denominator;
  
  double sum1 = -sum(Z % log(p_bag) + (1 - Z) % log(1 - p_bag));
  return sum1;
}
