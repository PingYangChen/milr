#ifndef common_h
#define common_h

#include <RcppArmadillo.h>
using namespace arma;

void chk_mat(const mat& x, const std::string& varName);
arma::vec logit(const arma::mat& X, const arma::vec& beta);
#endif
