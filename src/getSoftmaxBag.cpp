#include "common.h"

//' Get bag response function via softmax approach
//'
//' Get the class of bags via softmax approach.
//' 
//' @param X A matrix, the design matrix.
//' @param beta A vector, the coefficients.
//' @param bag A vector, the id of bags.
//' @return A vector. The classes of bags.
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector getSoftmaxBag(const arma::mat& X, const arma::vec& beta, const arma::vec& bag, 
                                  const double& alpha){
  chk_mat(X, "X");
  chk_mat(beta, "beta");
  chk_mat(bag, "bag");
  
  arma::uvec bag2 = arma::conv_to<arma::uvec>::from(bag - 1);
  arma::uvec uniBag = sort(arma::unique(bag2)), idx;
  arma::vec p1, tmp;
  double prob = 0.0;
  Rcpp::IntegerVector out(uniBag.n_elem);
  
  for (arma::uword i = 0; i < uniBag.n_elem; ++i) {
    idx = arma::find(bag2 == uniBag(i));
    p1 = logit(X.rows(idx), beta);
    tmp = exp(p1.elem(arma::find_finite(p1)) * alpha);
    prob = sum(p1.elem(arma::find_finite(p1)) % tmp) / sum(tmp);
    out[i] = prob > 0.5 ? 1 : 0;
  }
  return out;
}
