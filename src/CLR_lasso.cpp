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
arma::vec logit(const arma::mat& X, const arma::vec& beta){
  return pow(1.0 + exp(-X*beta), -1.0);
}

// q is the expected value of the instance (Y) given that the label (Z) is 1
arma::vec q(const arma::vec& p_Y, const arma::uvec& ID){
  // p_bag is P(label of bag_i is 1) = 1-prod(1-p_ij)
  uvec ID_table = unique(ID); 
  vec p_bag = ones<vec>(ID_table.n_elem);
  // assuming that the bags(ID) are labeled from 1 to the number of bags and are sorted
  for(uword i = 0; i < ID_table.size(); i++)
    p_bag(ID_table(i)) = prod(1.0 - p_Y.elem(find(ID == ID_table(i))));
  p_bag = 1 - p_bag;
  
  vec q = p_Y;
  q.elem(find(p_Y == 0.0)).zeros();
  uvec loc_nonzero_q = find(q > 0.0);
  q.elem(loc_nonzero_q) /= p_bag.elem(ID.elem(loc_nonzero_q));
  q.elem(find_nonfinite(q)).ones();
  return(q);
}

//[[Rcpp::export]]
arma::vec CLR_lasso(const arma::vec& Z, const arma::mat& X, const arma::vec& ID_dbl, 
                    const arma::vec& init_beta, const double& lambda, double alpha = 1,
                    double maxit = 500){
  double iter = 1.0;
  uword p = X.n_cols, n = X.n_rows;
  uvec ID = conv_to<uvec>::from(ID_dbl);
  ID -= 1;
  double diff = 1.0, W = 0.25;
  // use the upper bound 0.25 to approximate W = p(1-p)
  vec beta = init_beta;
  
  // X is normalized prior to data analysis
  // so XWX = W * sum(x^2) = W * (n-1)
  double XWX = W * ((double) n - 1.0);
  vec new_beta(p), p_vec(n), q_vec(n), U(n);
  
  while(diff > 1e-5 && iter < maxit){
    p_vec = logit(X, beta);
    // To avoid coefficients diverging in order to achieve fitted probabilities of 0 
    // or 1, when a probability is within 10^(-5) of 1, we set it to 1. 0 is treated
    // similarly.
    p_vec.elem(find(p_vec < 1e-5)).zeros();
    q_vec = q(p_vec, ID);
    vec S = trans(X) * (Z % q_vec - p_vec);
    for (uword k = 0; k < p; k++)
    {
      double tmp = S[k] + XWX * beta[k];
      if (k == 0)
      {
        new_beta[k] = tmp/XWX;
      } else
      {
        if(abs(tmp) <= lambda*alpha)
          new_beta[k] = 0;
        if(tmp > lambda*alpha)
          new_beta[k] = (tmp-lambda)/(XWX+lambda*(1-alpha));
        if(tmp < -lambda*alpha)
          new_beta[k] = (tmp+lambda)/(XWX+lambda*(1-alpha));
      }
    }
    // if the relative difference is less than tol, stop iterating
    diff = norm(new_beta - beta, 2) / norm(beta, 2);
    beta = new_beta;
    iter++;
  }
  return(beta);
}

