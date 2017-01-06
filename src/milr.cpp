#include "common.h"

//[[Rcpp::export]]
arma::vec EM_Y(const arma::field<arma::uvec>& bagField, const arma::vec& p_instance) {
  // p_bag is P(label of bag_i is 1) = 1-prod(1-p_ij)
  vec p_bag = ones<vec>(p_instance.n_elem);
  for (uword i =0; i < bagField.n_elem; ++i) 
    p_bag.elem(bagField(i)).fill(1 - prod(1.0 - p_instance.elem(bagField(i))));
  
  vec q = p_instance;
  q.elem(find(p_instance == 0.0)).zeros();
  uvec loc_nonzero_q = find(q > 0.0);
  q.elem(loc_nonzero_q) /= p_bag.elem(loc_nonzero_q);
  q.elem(find_nonfinite(q)).ones();
  return(q);
}

//[[Rcpp::export]]
arma::vec milr_cpp(const arma::vec& Z, const arma::mat& X, const arma::vec& bag,
                   const arma::vec& init_beta, const double& lambda, const double& tolerance, 
                   const double& alpha,  const double& maxit) {
  chk_mat(Z, "Z");
  chk_mat(X, "X");
  chk_mat(bag, "bag");
  chk_mat(init_beta, "init_beta");
  
  uword p = X.n_cols, n = X.n_rows;
  // convert bag to uword vec
  uvec bag2 = conv_to<uvec>::from(bag);
  uvec uniBag = sort(unique(bag2));
  arma::field<arma::uvec> bagField(uniBag.n_elem);
  for (uword i = 0; i < uniBag.n_elem; ++i)
    bagField(i) = find(bag2 == uniBag(i));
  
  // use the upper bound 0.25 to approximate W = p(1-p)
  double iter = 1.0, eps = 1.0, W = 0.25;
  
  // X is normalized prior to data analysis
  // so XWX = W * sum(x^2) = W * (n-1)
  double XWX = W * ((double) n - 1.0);
  vec beta = init_beta, new_beta(p), p_vec(n), q_vec(n), S(p);
  
  while (eps > tolerance && iter < maxit) {
    p_vec = logit(X, beta);
    // To avoid coefficients diverging in order to achieve fitted probabilities of 0
    // or 1, when a probability is within 10^(-5) of 1, we set it to 1. 0 is treated
    // similarly.
    p_vec.elem(find(p_vec < 1e-5)).zeros();
    q_vec = EM_Y(bagField, p_vec);
    S = X.t() * (Z % q_vec - p_vec);
    for (uword k = 0; k < p; ++k) {
      double tmp = S[k] + XWX * beta[k];
      if (k == 0) {
        new_beta[k] = tmp / XWX;
      } else {
        if(abs(tmp) <= lambda * alpha)
          new_beta[k] = 0.0;
        if(tmp > lambda * alpha)
          new_beta[k] = (tmp - lambda) / (XWX + lambda * (1.0 - alpha));
        if(tmp < -lambda * alpha)
          new_beta[k] = (tmp + lambda) / (XWX + lambda * (1.0 - alpha));
      }
    }
    // if the relative difference is less than tol, stop iterating
    eps = norm(new_beta - beta, 2) / norm(beta, 2);
    beta = new_beta;
    ++iter;
  }
  return beta;
}
