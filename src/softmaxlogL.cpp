#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double softmaxlogL(arma::vec bag, arma::mat instance, arma::vec label_bag, arma::vec beta, double alpha){
  vec p_instance = pow(1.0 + exp(-X*beta), -1.0);
  p_instance.elem(find(p_instance < 1e-7)).fill(1e-7);
  p_instance.elem(find(p_instance > 1-1e-7)).fill(1-1e-7);
  vec p_bag = zeros<vec>(label_bag.n_elem);
  vec denominator = zeros<vec>(label_bag.n_elem);
  
  for(uword i = 0; i < bag.n_elem; i++)
  {
    p_bag[bag[i] - 1] += p_instance[i] * exp(alpha * p_instance[i]);
    denominator[bag[i] - 1] += exp(alpha * p_instance[i]);
  }
  p_bag /= denominator;
  
  double sum1 = -sum(label_bag % log(p_bag) + (1 - label_bag) % log(1-p_bag));
  return sum1;
}