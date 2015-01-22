#ifndef COMP
#define COMP
using namespace Rcpp;

extern double com_compute_log_z(double lambda, double nu, double log_error, int maxit);
extern double com_compute_z(double lambda, double nu, double log_error, int maxit);
extern double logsumexp(NumericVector x);

extern NumericVector dcom(NumericVector x, double lambda, double nu, double z, bool log, double log_error, int maxit);
extern NumericVector dcom_parallel(NumericVector x, double lambda, double nu, double z, bool log, double, int maxit);
extern double dcom_single(double x, double lambda, double nu, double z);

extern NumericVector pcom(NumericVector q, double lambda, double nu, double z, bool log, double log_error, int maxit);
extern NumericVector pcom_parallel(NumericVector q, double lambda, double nu, double z, bool log, double, int maxit);
extern double pcom_single(double q, double lambda, double nu, double z);

extern NumericVector qcom(NumericVector p, double lambda, double nu, double z, bool log, double log_error, int maxit);
extern NumericVector qcom_parallel(NumericVector p, double lambda, double nu, double z, bool log, double, int maxit);
extern double qcom_single(double p, double lambda, double nu, double z);

extern NumericVector rcom(int n, double lambda, double nu, double z, bool log, double log_error, int maxit);
#endif
