#include <Rcpp.h>
#include <math.h>
#include <RcppParallel.h>
#include <boost/bind.hpp>
#include "compoisson.h"
#include "parallel-workers.h"

using namespace RcppParallel;


// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel, BH)]]

//' @import RcppParallel
//' @export
// [[Rcpp::export]]
NumericVector dcom_parallel(NumericVector x, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {
  
  // Perform argument checking
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  double log_z;

  if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error, maxit);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector d(x.size());
  
  Dcom ddcom(x, lambda, nu, log_z, d);
  
  parallelFor(0, x.size(), ddcom);
  
  if(!log) {
    d = Rcpp::exp(d);
  }

  return d;
}

//' @import RcppParallel
//' @export
// [[Rcpp::export]]
NumericVector pcom_parallel(NumericVector q, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {
  
  // Perform argument checking
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  double log_z;

  if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error, maxit);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector p(q.size());
  
  Pcom ppcom(q, lambda, nu, log_z, p);
  
  parallelFor(0, q.size(), ppcom);
  
  if(log) {
    p = Rcpp::log(p);
  }

  return p;
}




double pcom_single(double q, double lambda, double nu, double log_z) {
  double p = 0;
    for (int i = 0; i <= q; ++i) {
      p += exp(dcom_single(i, lambda, nu, log_z));
    }
  return p;
  }


