#include <RcppArmadillo.h>
#include <math.h>
#include "cmp.h"
#include "parallel-workers.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo]]
// [[Rcpp::interfaces(r, cpp)]]

//' @export
// [[Rcpp::export]]
double cmp_loglik(NumericMatrix x, NumericVector lambda, NumericVector nu, NumericVector z = NumericVector::create(NA), double log_error_z = 1e-6, int maxit_z = 10000) { //, bool parallel = false) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec lls = Rcpp::as<arma::vec>(dcmp(x(_, 0), lambda, nu, z, true, log_error_z, maxit_z)); //, parallel));
  return arma::dot(counts2, lls);
}



// //' @export
// // [[Rcpp::export]]
// IntegerVector cmp_loglik2(NumericVector x, double lambda, double nu, double z = NA_REAL, double log_error_z = 1e-6, int maxit_z = 10000, bool parallel = false) {
//   return(table(x));  
// }

//'@export
// [[Rcpp::export]]
double pois_loglik(NumericMatrix x, double lambda) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec lls(counts.size());
  for(int i = 0; i < lls.size(); ++i) {
    lls[i]= R::dpois(x(i,0), lambda, 1);
  }
  return arma::dot(counts2, lls);
}

//'@export
// [[Rcpp::export]]
double nb_loglik(NumericMatrix x, double mu, double size) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec lls(counts.size());
  for(int i = 0; i < lls.size(); ++i) {
    lls[i]= Rf_dnbinom_mu(x(i,0), size, mu, 1);
  }  return arma::dot(counts2, lls);
}

//'@export
// [[Rcpp::export]]
double pb_loglik(NumericMatrix x, double size, double prob) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec lls(counts.size());
  for(int i = 0; i < lls.size(); ++i) {
    lls[i]= Rf_dbinom(x(i,0), size, prob, 1);
  }  return arma::dot(counts2, lls);
}


//' @export
// [[Rcpp::export]]
double cmp_kld(NumericMatrix x, NumericVector lambda, NumericVector nu, NumericVector z = NumericVector::create(NA), double log_error_z = 1e-6, int maxit_z = 10000) { //, bool parallel = false) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec counts3 = counts2;
  counts3.elem(arma::find(counts2 <=0)).zeros();
  counts2.elem(arma::find(counts2 <=0)).ones();
  arma::vec lls = Rcpp::as<arma::vec>(dcmp(x(_, 0), lambda, nu, z, true, log_error_z, maxit_z)); //, parallel));
  return arma::dot(counts3, arma::log(counts2) - lls);
}

//'@export
// [[Rcpp::export]]
double pois_kld(NumericMatrix x, double lambda) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec counts3 = counts2;
  counts3.elem(arma::find(counts2 <= 0)).zeros();
  counts2.elem(arma::find(counts2 <= 0)).ones();
  arma::vec lls(counts.size());
  for(int i = 0; i < lls.size(); ++i) {
    lls[i]= R::dpois(x(i,0), lambda, 1);
  }
  return arma::dot(counts3, arma::log(counts2) - lls);
}

//'@export
// [[Rcpp::export]]
double nb_kld(NumericMatrix x, double mu, double size) {
  NumericVector counts = x(_,1);
  arma::vec counts2 = Rcpp::as<arma::vec>(counts);
  arma::vec counts3 = counts2;
  counts3.elem(arma::find(counts2 <= 0)).zeros();
  counts2.elem(arma::find(counts2 <= 0)).ones();
  arma::vec lls(counts.size());
  for(int i = 0; i < lls.size(); ++i) {
    lls[i]= Rf_dnbinom_mu(x(i,0), size, mu, 1);
  }
  return arma::dot(counts3, arma::log(counts2) - lls);
}
