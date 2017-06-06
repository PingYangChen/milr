#ifndef common_h
#define common_h

#include <RcppArmadillo.h>

void chk_mat(const arma::mat& x, const std::string& varName);
arma::vec logit(const arma::mat& X, const arma::vec& beta);
#endif
