#include <RcppArmadillo.h>
#include <math.h>
#include "compoisson.h"
#include "parallel-workers.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo]]
// [[Rcpp::interfaces(r, cpp)]]

//' @export
// [[Rcpp::export]]
double com_loglik(NumericMatrix x, double lambda, double nu, double z = NA_REAL, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec lls = Rcpp::as<arma::vec>(dcom(x(_, 0), lambda, nu, z, true, log_error_z, maxit_z, parallel));
  return arma::dot(counts2, lls);
}
