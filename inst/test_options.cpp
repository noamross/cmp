#include <Rcpp.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
using namespace Rcpp;

#define USE_RINTERNALS
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
int test_getoption(const char* opt) {
  //  const char* optname = as<char>(opt);
    const SEXP option_name = Rf_install(opt);
    return Rf_asInteger(Rf_GetOption1(option_name));
}
