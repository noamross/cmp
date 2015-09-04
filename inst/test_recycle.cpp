#include <Rcpp.h>
#include <math.h>
#include "cmp.h"
#include "parallel-workers.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
NumericVector dcmp(NumericVector x, NumericVector lambda, NumericVector nu, NumericVector z = NumericVector::create(NA_REAL), bool log = false, double log_error_z = 1e-6, int maxit_z = 10000) { //, bool parallel = false) {
  
  
  int size = max(NumericVector::create(x.size(), lambda.size(), nu.size(), z.size()));
  
  x = rep_len(x, size);
  lambda = rep_len(lambda, size);
  nu = rep_len(nu, size);
  z = rep_len(z, size);
  
  NumericVector d(size);
  NumericVector log_z(size);
  
  for (int i = 0; i < size; i++) {
    if(ISNA(z[i])) {
      log_z[i] = compute_log_z(lambda[i], nu[i], log_error_z, maxit_z);
    } else{
      log_z[i] = std::log(z[i]);
    }
  }
  
//   if(parallel) {
//     Dcmp ddcmp(x, lambda, nu, log_z, d);
//     parallelFor(0, x.size(), ddcmp);
//     
//   } else {
    for (int i = 0; i < x.size(); ++i) {  
      d[i] = dcmp_single(x[i], lambda[i], nu[i], log_z[i]);
    }
  // }
  
  if(!log) {
    d = Rcpp::exp(d);
  }
  
  return d;
}


double dcmp_single(double x, double lambda, double nu, double log_z) {
  double d;
  if (lambda < 0 || nu < 0) {
    d = NAN;
    if (x < 0 || x != floor(x)) {
      d = R_NegInf; 
    } else if (x == 0 && lambda == 0) {
      d = 0;
    } else {
      d = x * std::log(lambda) - nu * Rcpp::internal::lfactorial(x) - log_z;
    }
    return d;
  }