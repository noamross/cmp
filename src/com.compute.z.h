#ifndef COMPZ
#define COMPZ
using namespace Rcpp;

extern double compute_log_z(double lambda, double nu, double log_error, int maxit);
extern double logsumexp(NumericVector x);
extern double llf(NumericVector x);
#endif
