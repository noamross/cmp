#ifndef COMP
#define COMP
using namespace Rcpp;

extern double com_compute_log_z(double lambda, double nu, double log_error, int maxit);
extern double com_compute_z(double lambda, double nu, double log_error, int maxit);
extern double logsumexp(NumericVector x);
extern double dcom(double x, double lambda, double nu, double z, double log_error, int maxit);
extern double pcom(double q, double lambda, double nu, bool log, double z, double log_error, int maxit);
extern int qcom(double p, double lambda, double nu, bool log, double z, double log_error, int maxit);
extern NumericVector rcom(int n, double lambda, double nu, bool log, double z, double log_error, int maxit);
#endif
