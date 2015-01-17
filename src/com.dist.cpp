#include <Rcpp.h>
#include "compoisson.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' @export
// [[Rcpp::export]]
double d_com(double x, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
    // Perform argument checking
  if (lambda < 0 || nu < 0) {
		Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	} else if (x < 0 || x != floor(x)) {
		return (0);
	} else if (ISNA(z) && !log) {
    z = compute_z(lambda, nu, log_error, maxit);
	} else if (ISNA(z) && log) {
    z = compute_log_z(lambda, nu, log_error, maxit);
	}  
  if (log) {
   return ((x * std::log(lambda) - nu * (Rcpp::internal::lfactorial(x)) - z));
  } else {
    return (pow(lambda, x) * (pow(Rcpp::internal::factorial(x), -nu) / z));
  }
}

//' @export
// [[Rcpp::export]]
double p_com(double q, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  double prob = 0;
  for (int i = 0; i <= q; ++i) {
    prob += d_com(i, lambda, nu, false, z, log_error, maxit);
  }
  if(log) {
    prob = std::log(prob);
  }
  return(prob);
  }

//' @export
// [[Rcpp::export]]
int q_com(double p, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  
  double prob = 0;
  int i = 0;
  if (log) {
    p = exp(p);
  }
  while (prob < p) {
    prob += d_com(i, lambda, nu, false, z, log_error, maxit);
    i += 1;
  }
  
  return (i - 1);
  }

//' @export
// [[Rcpp::export]]
NumericVector r_com(int n, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
    NumericVector vals(n);
    for(int i = 0; i < n; ++i) {
      vals[i] = q_com(as<double>(runif(1)), lambda, nu, log, z, log_error, maxit);
    }
    return vals;
}
