#include <Rcpp.h>
#include "compoisson.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' The COM-Poisson Distribution
//'
//' Probability mass function and random generation for the COM-Poisson
//' distribution for given values of the parameters.
//'
//' Computes the probability mass function of the COM-Poisson distribution
//' \deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}\deqn{ f(x) =
//' \frac{1}{Z(\lambda,\nu)}\frac{\lambda^x}{(x!)^\nu}. }{f(x) = (1/Z)
//' (lambda^x)/(x!^nu).}\deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}
//'
//' @aliases dcom pcom qcom rcom 
//' @param x level to evaluate the PMF at
//' @param lambda value of lambda parameter
//' @param nu value of nu parameter
//' @param z normalizing constant, computed if not specified
//' @param n number of random values to return
//' @param log.z natural log of z
//' @return \code{dcom} gives the probability that a random COM-Poisson variable
//' X takes value x.  \code{rcom} gives a vector of \code{n} random values
//' sampled from the COM-Poisson distribution.
//' @author Jeffrey Dunn
//' @seealso \code{\link{com.loglikelihood}}, \code{\link{com.log.density}}
//' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
//' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
//' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
//' v54, pp. 127-142, 2005.
//' @keywords models
//' @export
// [[Rcpp::export]]
double dcom(double x, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  
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
double pcom(double q, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  double prob = 0;
  for (int i = 0; i <= q; ++i) {
    
    prob += dcom(i, lambda, nu, false, z, log_error, maxit);
  }
  if(log) {
    
    prob = std::log(prob);
  }
  return(prob);
}

//' @export
// [[Rcpp::export]]
int qcom(double p, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  
  double prob = 0;
  int i = 0;
  if (log) {
    
    p = exp(p);
  }
  while (prob < p) {
    prob += dcom(i, lambda, nu, false, z, log_error, maxit);
    i += 1;
  }
  
  return (i - 1);
}

//' @export
// [[Rcpp::export]]
NumericVector rcom(int n, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit=100) {
  
  NumericVector vals(n);
  for(int i = 0; i < n; ++i) {
    vals[i] = qcom(as<double>(runif(1)), lambda, nu, log, z, log_error, maxit);
  }
  
  return vals;
}
