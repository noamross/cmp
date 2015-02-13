#include <Rcpp.h>
#include <math.h>
#include "cmp.h"
#include "parallel-workers.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::interfaces(r, cpp)]]

double dcmp_single(double x, double lambda, double nu, double log_z) {
  double d;
  if (x < 0 || x != floor(x)) {
    d = R_NegInf; 
  } else {
    d = x * std::log(lambda) - nu * Rcpp::internal::lfactorial(x) - log_z;
  }
  return d;
  }

double pcmp_single(double q, double lambda, double nu, double log_z) {
  double p = 0;
    for (int i = 0; i <= q; ++i) {
      p = logsumexp(p, dcmp_single(i, lambda, nu, log_z));
    }
  return p;
  }

double qcmp_single(double p, double lambda, double nu, double log_z) {
  double q;

  if (p == R_NegInf) {
    q = 0;
  } else if (p==0) {
    q = R_PosInf;
  } else {
    double prob = -std::numeric_limits<double>::infinity();
    q = 0;
    while (prob < p) {
     // Rcout << prob << std::endl;
      prob = logsumexp(prob, dcmp_single(q, lambda, nu, log_z));
      q += 1;
    }
    q = q - 1;
  }
  return q;
}




//' The Conway-Maxwell-Poisson Distribution
//'
//' Probability mass function and random generation for the Conway-Maxwell-
//' Poisson Distribution distribution for given values of the parameters.
//'
//' Computes the probability mass function of the CMP distribution
//' \deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}\deqn{ f(x) =
//' \frac{1}{Z(\lambda,\nu)}\frac{\lambda^x}{(x!)^\nu}. }{f(x) = (1/Z)
//' (lambda^x)/(x!^nu).}\deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}
//'
//' @aliases dcmp pcmp qcmp rcmp 
//' @param x level to evaluate the PMF at
//' @param lambda value of lambda parameter
//' @param nu value of nu parameter
//' @param z normalizing constant, computed if not specified
//' @param n number of random values to return
//' @param log.z natural log of z
//' @return \code{dcmp} gives the probability that a random CMP variable
//' X takes value x.  \code{rcmp} gives a vector of \code{n} random values
//' sampled from the CMP distribution.
//' @author Jeffrey Dunn
//' @seealso \code{\link{cmp_loglikelihood}},
//' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
//' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
//' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
//' v54, pp. 127-142, 2005.
//' @keywords models
//' @export
// [[Rcpp::export]]
NumericVector dcmp(NumericVector x, double lambda, double nu, double z = NA_REAL, bool log_p = false, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
  
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  double log_z;
  
  if (ISNA(z)) {
    log_z = compute_log_z(lambda, nu, log_error_z, maxit_z);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector d(x.size());
  
  if(parallel) {
    Dcmp ddcmp(x, lambda, nu, log_z, d);
    parallelFor(0, x.size(), ddcmp);
    
  } else {
    for (int i = 0; i < x.size(); ++i) {  
      d[i] = dcmp_single(x[i], lambda, nu, log_z);
    }
  }
  
  if(!log_p) {
    d = Rcpp::exp(d);
  }
  
  return d;
}


//' @rdname dcmp
//' @export
// [[Rcpp::export]]
NumericVector pcmp(NumericVector q, double lambda, double nu, double z = NA_REAL, bool log_p = false, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
  
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
  
  double log_z;
  
  if (ISNA(z)) {
    log_z = compute_log_z(lambda, nu, log_error_z, maxit_z);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector p(q.size());
  
  if(parallel) {
    Pcmp ppcmp(q, lambda, nu, log_z, p);
    parallelFor(0, q.size(), ppcmp);  
  } else {
    for (int i = 0; i < q.size(); ++i) {
      p[i] = pcmp_single(q[i], lambda, nu, log_z);
    }
  }    
  if(!log_p) {
    p = Rcpp::exp(p);
  }
  return(p);
}


//' @rdname dcmp
//' @export
// [[Rcpp::export]]
NumericVector qcmp(NumericVector p, double lambda, double nu, double z = NA_REAL, bool log_p = false, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
  
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  } 
    
  double log_z;
  
  if (ISNA(z)) {
    log_z = compute_log_z(lambda, nu, log_error_z, maxit_z);
  } else {
    log_z = std::log(z);
  }
  
  NumericVector lp;
  if (!log_p) {
    lp = log(p);
  } else {
    lp = p;
  }

  NumericVector q(lp.size());
  
  if(parallel) {  
    Qcmp qqcmp(p, lambda, nu, log_z, q);
    parallelFor(0, p.size(), qqcmp);
  } else {
    for (int i = 0; i < q.size(); ++i) {
      q[i] = qcmp_single(lp[i], lambda, nu, log_z);
    }
  }
  
  return q;
}

//' @rdname dcmp
//' @export
// [[Rcpp::export]]
NumericVector rcmp(int n, double lambda, double nu, double z = NA_REAL, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
    return(qcmp(runif(n, 0, 1), lambda, nu, z, false, log_error_z, maxit_z, parallel));
}
