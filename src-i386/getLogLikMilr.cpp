#include "common.h"
#include <RcppParallel.h>

struct Worker_getLogLikMilr : public RcppParallel::Worker {
  const arma::uvec& bag2;
  const arma::uvec& uniBag;
  const arma::vec& y;
  const arma::mat& X;
  const arma::vec& beta;
  double logLikMilr;
  
  Worker_getLogLikMilr(const arma::uvec& bag2, const arma::uvec& uniBag, const arma::vec& y, 
                       const arma::mat& X, const arma::vec& beta)
    : bag2(bag2), uniBag(uniBag), y(y), X(X), beta(beta), logLikMilr(0) {}
  
  Worker_getLogLikMilr(const Worker_getLogLikMilr& getLogLikMilr_worker, RcppParallel::Split)
    : bag2(getLogLikMilr_worker.bag2), uniBag(getLogLikMilr_worker.uniBag), y(getLogLikMilr_worker.y), 
    X(getLogLikMilr_worker.X),  beta(getLogLikMilr_worker.beta), logLikMilr(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (arma::uword i = begin; i < end; ++i) {
      arma::uvec idx = arma::find(bag2 == uniBag(i));
      double prob = std::min(1 - 1e-16, std::max(1e-16, 1 - prod(1.0 - logit(X.rows(idx), beta))));
      logLikMilr += y(idx(0)) * std::log(prob) + (1 - y(idx(0))) * std::log(1 - prob);
    }
  }
  
  // join my value with that of another Sum
  void join(const Worker_getLogLikMilr& rhs) {
    logLikMilr += rhs.logLikMilr; 
  }
};

// [[Rcpp::export]]
double getLogLikMilr(const arma::vec& beta, const arma::vec& y, 
                     const arma::mat& X, const arma::vec& bag) {
  chk_mat(beta, "beta");
  chk_mat(y, "y");
  chk_mat(X, "X");
  chk_mat(bag, "bag");
  
  arma::uvec bag2 = arma::conv_to<arma::uvec>::from(bag);
  arma::uvec uniBag = arma::sort(arma::unique(bag2)), idx;
  Worker_getLogLikMilr getLogLikMilr_worker(bag2, uniBag, y, X, beta);
  parallelReduce(0, uniBag.n_elem, getLogLikMilr_worker);
  return getLogLikMilr_worker.logLikMilr;
}
