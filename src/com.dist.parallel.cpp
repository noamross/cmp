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

