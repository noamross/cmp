#include <Rcpp.h>
#include "compoisson.h"


using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double compute_z(double lambda, double nu, double log_error = 0.001, int maxit=100) {
  return(exp(compute_log_z(lambda, nu, log_error, maxit)));
}

//' @export
// [[Rcpp::export]]
double compute_log_z(double lambda, double nu, double log_error = 0.001, int maxit=100)
{
  // Perform argument checking
	if (lambda < 0 || nu < 0) {
		Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	}
	// Initialize values
	double z = R_NegInf;
	double z_last = 0;
	int j = 0;
	while ((std::abs(z - z_last) > log_error) && j <= maxit)
	{
		z_last = z;
		z = logsumexp(NumericVector::create(z, j * log(lambda) - nu * Rcpp::internal::lfactorial(j)));
		j += 1;
	}
	return z;
}

//' @export
// [[Rcpp::export]]
double logsumexp(NumericVector x) {
    return(log(sum(exp(x - max(x)))) + max(x));
}
