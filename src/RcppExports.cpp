// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/compoisson.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// com_compute_z
double com_compute_z(double lambda, double nu, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_com_compute_z_try(SEXP lambdaSEXP, SEXP nuSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        double __result = com_compute_z(lambda, nu, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_com_compute_z(SEXP lambdaSEXP, SEXP nuSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_com_compute_z_try(lambdaSEXP, nuSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// com_compute_log_z
double com_compute_log_z(double lambda, double nu, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_com_compute_log_z_try(SEXP lambdaSEXP, SEXP nuSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        double __result = com_compute_log_z(lambda, nu, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_com_compute_log_z(SEXP lambdaSEXP, SEXP nuSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_com_compute_log_z_try(lambdaSEXP, nuSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// logsumexp
double logsumexp(NumericVector x);
static SEXP compoisson_logsumexp_try(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        double __result = logsumexp(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_logsumexp(SEXP xSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_logsumexp_try(xSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// dcom
NumericVector dcom(NumericVector x, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_dcom_try(SEXP xSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< bool >::type log(logSEXP );
        Rcpp::traits::input_parameter< double >::type z(zSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        NumericVector __result = dcom(x, lambda, nu, log, z, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_dcom(SEXP xSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_dcom_try(xSEXP, lambdaSEXP, nuSEXP, logSEXP, zSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// pcom
NumericVector pcom(NumericVector q, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_pcom_try(SEXP qSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP );
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< bool >::type log(logSEXP );
        Rcpp::traits::input_parameter< double >::type z(zSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        NumericVector __result = pcom(q, lambda, nu, log, z, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_pcom(SEXP qSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_pcom_try(qSEXP, lambdaSEXP, nuSEXP, logSEXP, zSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// qcom
NumericVector qcom(NumericVector p, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_qcom_try(SEXP pSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< bool >::type log(logSEXP );
        Rcpp::traits::input_parameter< double >::type z(zSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        NumericVector __result = qcom(p, lambda, nu, log, z, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_qcom(SEXP pSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_qcom_try(pSEXP, lambdaSEXP, nuSEXP, logSEXP, zSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// rcom
NumericVector rcom(int n, double lambda, double nu, bool log = false, double z = NA_REAL, double log_error = 0.001, int maxit = 100);
static SEXP compoisson_rcom_try(SEXP nSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< bool >::type log(logSEXP );
        Rcpp::traits::input_parameter< double >::type z(zSEXP );
        Rcpp::traits::input_parameter< double >::type log_error(log_errorSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        NumericVector __result = rcom(n, lambda, nu, log, z, log_error, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP compoisson_rcom(SEXP nSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP logSEXP, SEXP zSEXP, SEXP log_errorSEXP, SEXP maxitSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(compoisson_rcom_try(nSEXP, lambdaSEXP, nuSEXP, logSEXP, zSEXP, log_errorSEXP, maxitSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}

// validate (ensure exported C++ functions exist before calling them)
static int compoisson_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*com_compute_z)(double,double,double,int)");
        signatures.insert("double(*com_compute_log_z)(double,double,double,int)");
        signatures.insert("double(*logsumexp)(NumericVector)");
        signatures.insert("NumericVector(*dcom)(NumericVector,double,double,bool,double,double,int)");
        signatures.insert("NumericVector(*pcom)(NumericVector,double,double,bool,double,double,int)");
        signatures.insert("NumericVector(*qcom)(NumericVector,double,double,bool,double,double,int)");
        signatures.insert("NumericVector(*rcom)(int,double,double,bool,double,double,int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP compoisson_RcppExport_registerCCallable() { 
    R_RegisterCCallable("compoisson", "compoisson_com_compute_z", (DL_FUNC)compoisson_com_compute_z_try);
    R_RegisterCCallable("compoisson", "compoisson_com_compute_log_z", (DL_FUNC)compoisson_com_compute_log_z_try);
    R_RegisterCCallable("compoisson", "compoisson_logsumexp", (DL_FUNC)compoisson_logsumexp_try);
    R_RegisterCCallable("compoisson", "compoisson_dcom", (DL_FUNC)compoisson_dcom_try);
    R_RegisterCCallable("compoisson", "compoisson_pcom", (DL_FUNC)compoisson_pcom_try);
    R_RegisterCCallable("compoisson", "compoisson_qcom", (DL_FUNC)compoisson_qcom_try);
    R_RegisterCCallable("compoisson", "compoisson_rcom", (DL_FUNC)compoisson_rcom_try);
    R_RegisterCCallable("compoisson", "compoisson_RcppExport_validate", (DL_FUNC)compoisson_RcppExport_validate);
    return R_NilValue;
}
