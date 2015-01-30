#ifndef COMP
#define COMP

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;

extern double com_compute_z(double lambda, double nu, double log_error_z, int maxit_z);
extern double com_compute_log_z(double lambda, double nu, double log_error_z, int maxit_z);
extern double com_compute_log_z_approx(double lambda, double nu);
extern double logsumexp(NumericVector x);
extern double logsumexp(double x, double y);
extern double logdiffexp(double x, double y);

extern NumericVector dcom(NumericVector x, double lambda, double nu, double z, bool log_p, double log_error_z, int maxit_z, bool parallel);
extern double dcom_single(double x, double lambda, double nu, double z);

extern NumericVector pcom(NumericVector q, double lambda, double nu, double z, bool log_p, double log_error_z, int maxit_z, bool parallel);
extern double pcom_single(double q, double lambda, double nu, double z);

extern NumericVector qcom(NumericVector p, double lambda, double nu, double z, bool log_p, double log_error_z, int maxit_z, bool parallel);
extern double qcom_single(double p, double lambda, double nu, double z);

extern NumericVector rcom(int n, double lambda, double nu, double z, double log_error_z, int maxit_z, bool parallel);

extern double com_mean(double lambda, double nu, double log_error, int maxit, double z, double log_error_z, int maxit_z, bool parallel);
extern double com_log_mean(double lambda, double nu, double log_error, int maxit, double z, double log_error_z, int maxit_z, bool parallel);
extern double com_mean_approx(double lambda, double nu);
extern double com_log_mean_approx(double lambda, double nu);

extern double com_var(double lambda, double nu, double log_error, int maxit, double z, double log_error_z, int maxit_z, bool parallel);
extern double com_log_var(double lambda, double nu, double log_error, int maxit, double z, double log_error_z, int maxit_z, bool parallel);
extern double com_var_approx(double lambda, double nu);
extern double com_log_var_approx(double lambda, double nu);

extern double com_loglik(NumericVector x, double lambda, double nu, double z, double log_error_z, int maxit_z, bool parallel);
extern double pois_loglik(NumericMatrix x, double lambda);
extern double nb_loglik(NumericMatrix x, double prob, double mu);  

extern double com_compute_log_z_old(double lambda, double nu, double log_error_z, int maxit_z);

#endif
