#include <Rcpp.h>
#include <math.h>
#include "compoisson.h"
#include "parallel-workers.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Computes Expectation of a Function of a COM-Poisson Random Variable
//'
//' Computes an expectation of a function of a COM-Poisson random variable.
//'
//' Computes the expectation \eqn{E[f(X)]}{E[f(X)]} where X is a COM-Poisson
//' random variable.
//'
//' @param f function taking as a single argument the value of x
//' @param lambda value of lambda parameter
//' @param nu value of nu parameter
//' @param log.error precision in the log of the expectation
//' @return The expectation as a real number.
//' @author Jeffrey Dunn
//' @seealso \code{\link{com.mean}}, \code{\link{com.var}},
//' \code{\link{com.fit}}
//' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
//' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
//' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
//' v54, pp. 127-142, 2005.
//' @keywords models
//' @export
// [[Rcpp::export]]
double com_mean(double lambda, double nu,
                       double log_error = 1e-6, int maxit=1e6,
                       double z = NA_REAL, double log_error_z = 1e-6,
                       int maxit_z = 10000, bool parallel = false) {
  
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  if(exp(log(lambda)/nu) > maxit_z) return(pow(lambda, 1/nu) - (nu - 1)/(2*nu));
    
  double log_z;
  
   if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error_z, maxit_z);
  } else {
    log_z = std::log(z);
  }

  int j = 1;
  double ex = log(j) + dcom_single(j, lambda, nu, log_z);
  double last_ex = 0;

  j += 1; 
  
  while ((std::abs(last_ex - ex) > log_error) && j <= maxit) {
    last_ex = ex;
    ex = logsumexp(ex, log(j) + dcom_single(j, lambda, nu, log_z));
    j += 1;
  }
  
  return exp(ex);
}

//' @export
// [[Rcpp::export]]
double com_var(double lambda, double nu,
                       double log_error = 1e-6, int maxit=1e6,
                       double z = NA_REAL, double log_error_z = 1e-6,
                       int maxit_z = 10000, bool parallel = false) {
  
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
      
  double log_z;
  
   if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error_z, maxit_z);
  } else {
    log_z = std::log(z);
  }

  int j = 1;
  double ex = log(j) + dcom_single(j, lambda, nu, log_z);
  double last_ex = 0;

  j += 1; 
  
  while ((std::abs(last_ex - ex) > log_error) && j <= maxit) {
    last_ex = ex;
    ex = logsumexp(ex, log(pow(j,2)) + dcom_single(j, lambda, nu, log_z));
    j += 1;
  }
  
  ex = exp(ex) - pow(com_mean(lambda, nu, log_error, maxit, exp(log_z), log_error_z, maxit_z, parallel), 2);
  
  return ex;
}