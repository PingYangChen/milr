#include "common.h"
#include <RcppParallel.h>

struct Worker_getMilrProb: public RcppParallel::Worker {
  const uvec& bag2;
  const uvec& uniBag;
  const mat& X;
  const vec& beta;
  vec& prob;

  Worker_getMilrProb(const uvec& bag2, const uvec& uniBag, const mat& X, const vec& beta, vec& prob):
    bag2(bag2), uniBag(uniBag), X(X), beta(beta), prob(prob) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (uword i = begin; i < end; ++i) {
      uvec idx = find(bag2 == uniBag(i));
      prob(i) = 1 - prod(1.0 - logit(X.rows(idx), beta));
    }
  }
};

// [[Rcpp::export]]
arma::vec getMilrProb(const arma::vec& beta, const arma::mat& X, const arma::vec& bag) {
  chk_mat(beta, "beta");
  chk_mat(X, "X");
  chk_mat(bag, "bag");
  
  uvec bag2 = conv_to<uvec>::from(bag - 1);
  uvec uniBag = sort(unique(bag2)), idx;
  vec prob(uniBag.n_elem);
  Worker_getMilrProb getMilrProb_worker(bag2, uniBag, X, beta, prob);
  RcppParallel::parallelFor(0, uniBag.n_elem, getMilrProb_worker);
  return prob;
}
