#ifndef COMP
#define COMP
using namespace Rcpp;

extern double compute_log_z(double lambda, double nu, double log_error, int maxit);
extern double compute_z(double lambda, double nu, double log_error, int maxit);
extern double logsumexp(NumericVector x);
extern double d_com(double x, double lambda, double nu, double z, double log_error, int maxit);
extern double p_com(double q, double lambda, double nu, bool log, double z, double log_error, int maxit);
#endif
