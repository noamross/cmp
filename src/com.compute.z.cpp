#include <Rcpp.h>
#include "compoisson.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' Compute COM-Poisson Normalizing Constant
//'
//' Computes the normalizing constant in the COM-Poisson model for given values
//' of the parameters.
//'
//' \code{com_compute_z} computes the COM-Poisson normalizing constant \deqn{
//' }{z = Sum (lambda^j)/(j!^nu) }\deqn{ z = \sum_{i=0}^\infty
//' \frac{\lambda^j}{(j!)^\nu} }{z = Sum (lambda^j)/(j!^nu) }\deqn{ }{z = Sum
//' (lambda^j)/(j!^nu) } to the specified precision. If no precision is
//' specified, then the package default is used.
//'
//' \code{com_compute_log_z} is equivalent to \code{log(com_compute_z(lambda,
//' nu))} but provudes additional precision.
//'
//' @aliases com_compute_z com_compute_log_z
//' @param lambda Lambda value in COM-Poisson distribution
//' @param nu Nu value in COM-Poisson distribution
//' @param log.error Precision in the log of the normalizing constant
//' @return The normalizing constant as a real number with specified precision.
//' @author Jeffrey Dunn
//' @seealso \code{\link{com.fit}}
//' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
//' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
//' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
//' v54, pp. 127-142, 2005.
//' @keywords models
//' @examples
//'
//' data(insurance)
//' fit = com.fit(Lemaire)
//' z = com_compute_z(fit$lambda, fit$nu)
//'
//' @export
// [[Rcpp::export]]
double com_compute_z(double lambda, double nu, double log_error = 0.001, int maxit=1000) {
  return(exp(com_compute_log_z(lambda, nu, log_error, maxit)));
}


//' @rdname com_compute_z
//' @export
// [[Rcpp::export]]
double com_compute_log_z(double lambda, double nu, double log_error = 0.001, int maxit=1000) {
  // Perform argument checking
  if (lambda < 0 || nu < 0) {
    Rcpp::stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  }
  
  // Initialize values
  double z = R_NegInf;
  double z_last = 0;
  int j = 0;
  while ((std::abs(z - z_last) > log_error) && j <= maxit) {
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
