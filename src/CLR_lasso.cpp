#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

//' logit link function
//'
//' calculate the values of logit link
//' 
//' @param X A matrix, the design matrix.
//' @param beta A vector, the coefficients.
//' @return An vector of the values of logit link.
//' @examples
//' logit(matrix(c(1, 0.5, -0.2, 0.3), 2), c(-0.5, -3)) # 0.5249792, 0.2404891
//' @export
// [[Rcpp::export]]
arma::vec logit(arma::mat X, arma::vec beta){
  return pow(1 + exp(-X*beta), -1);
}

// q is the expected value of the instance (Y) given that the label (Z) is 1
arma::vec q(arma::vec p_Y, arma::vec ID){
  vec q = p_Y;
  // p_bag is P(label of bag_i is 1) = 1-prod(1-p_ij)
  vec ID_table = unique(ID); 
  vec p_bag = ones<arma::vec>(ID_table.n_elem);
  
  // assuming that the bags(ID) are labeled from 1 to the number of bags and are sorted
  for(uword i = 0; i < ID.size(); i++){
    p_bag[ID[i] - 1] *= 1-p_Y[i];
  }
  p_bag = 1-p_bag;
  
  for(uword j = 0; j < q.size(); j++){
    if(p_Y[j] == 0) q[j] = 0;
    else{
      double tmp = p_bag[ID[j]-1];
      if(tmp == 0){
        q[j] = 1;
      }
      else{
        q[j] = p_Y[j]/tmp;
      }
    }
  }
  return(q);
}

//[[Rcpp::export]]
arma::vec CLR_lasso(const arma::vec& Z, const arma::mat& X, const arma::vec& ID, 
                    const arma::vec& init_beta, double lambda, double alpha = 1,
                    double maxit = 500){
  double iter = 1.0;
  uword p = X.n_cols, n = X.n_rows;
  
  double tol = std::pow(0.1, 5.0), epsilon = std::pow(0.1, 5.0), 
    diff = 1.0, W = 0.25;
  // use the upper bound 0.25 to approximate W = p(1-p)
  vec beta = init_beta;
  // the initial value of beta is chosen uniformly at random in (0,1)
  
  vec new_beta(p), p_vec(n), q_vec(n), U(n); 
  // X is normalized prior to data analysis
  // so XWX = W * sum(x^2) = W * (n-1)
  double XWX = W * (n-1);
  
  while(diff > tol && iter < maxit){
    p_vec = logit(X, beta);
    // To avoid coefficients diverging in order to achieve fitted probabilities of 0 
    // or 1, when a probability is within 10^(-5) of 1, we set it to 1. 0 is treated
    // similarly.
    for(uword j = 0; j < p_vec.size(); j++){
      if(p_vec[j] < epsilon){
        p_vec[j] = 0;
      }
      else{
        if(p_vec[j] > 1-epsilon){
          p_vec[j] = 1;
        }
      }
    }
    q_vec = q(p_vec, ID);
    vec S = trans(X) * (Z % q_vec - p_vec);
    for(uword k = 0; k < p; k++){
      double tmp = S[k] + XWX * beta[k];
      if(k == 0){
        new_beta[k] = tmp/XWX;
      }
      else{
        if(abs(tmp) <= lambda){
          new_beta[k] = 0;
        }
        if(tmp > lambda){
          new_beta[k] = (tmp-lambda)/(XWX+lambda*(1-alpha));
        }
        if(tmp < -lambda){
          new_beta[k] = (tmp+lambda)/(XWX+lambda*(1-alpha));
        }
      }
    }
    // if the relative difference is less than tol, stop iterating
    diff = norm(new_beta - beta,2) / norm(beta, 2);
    beta = new_beta;
    iter++;
  }
  return(beta);
}

