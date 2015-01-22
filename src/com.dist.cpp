#include <Rcpp.h>
#include <math.h>
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
NumericVector dcom(NumericVector x, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {
  
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
  
  for (int i = 0; i < x.size(); ++i) {  
   d[i] = dcom_single(x[i], lambda, nu, log_z);
  }

  if(!log) {
    d = Rcpp::exp(d);
  }

  return d;
}

double dcom_single(double x, double lambda, double nu, double log_z) {
  double d;
  if (x < 0 || x != floor(x)) {
    d = R_NegInf; 
  } else {
    d = x * std::log(lambda) - nu * Rcpp::internal::lfactorial(x) - log_z;
  }
  return d;
  }

//' @rdname dcom
//' @export
// [[Rcpp::export]]
NumericVector pcom(NumericVector q, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {

  double log_z;

  // Perform argument checking
  // Perform argument checking
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } else if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error, maxit);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector p(q.size());
  double prob;
  
  for (int j = 0; j < q.size(); ++j) {
    prob = 0;
    for (int i = 0; i <= q[j]; ++i) {
      prob += exp(dcom_single(i, lambda, nu, log_z));
    }
    p[j] = prob;
  }
  
  if(log) {
    p = Rcpp::log(p);
  }
  return(p);
}

//' @rdname dcom
//' @export
// [[Rcpp::export]]
NumericVector qcom(NumericVector p, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {
    
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  double log_z;

  if (ISNA(z)) {
    log_z = com_compute_log_z(lambda, nu, log_error, maxit);
  } else {
    log_z = std::log(z);
  }

  if (log) {
    p = Rcpp::exp(p);
  }
  
  NumericVector q(p.size());
  double prob;
  int i;
  
  for (int j = 0; j < q.size(); ++j) {
    if (p[j] == 0) {
      q[j] = 0;
      continue;
    } else if (p[j] == 1) {
      q[j] = R_PosInf;
      continue;
    } else {
    prob = 0;
    i = 0;
    while (prob < p[j]) {
      prob += exp(dcom_single(i, lambda, nu, log_z));
      i += 1;
    }
    q[j] = i - 1;
    }
  }
  
  return q;
}

//' @rdname dcom
//' @export
// [[Rcpp::export]]
NumericVector rcom(int n, double lambda, double nu, double z = NA_REAL, bool log = false, double log_error = 0.001, int maxit=1000) {
    return(qcom(runif(n, 0, 1), lambda, nu, z, log, log_error, maxit));
}
