#ifndef PARWORK
#define PARWORK

#include <Rcpp.h>
#include <math.h>
#include <RcppParallel.h>
#include <boost/bind.hpp>
#include "cmp.h"

using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]


struct Dcmp : public Worker {   
   // source vector and parameters
   const RVector<double> input;
   const double lambda;
   const double nu;
   const double log_z;
   
   // destination vector
   RVector<double> output;
   
   // initialize with source and destination
   Dcmp(const NumericVector input, const double lambda, const double nu, const double log_z, NumericVector output) 
      : input(input), lambda(lambda), nu(nu), log_z(log_z), output(output) {}
   
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::boost::bind(&dcmp_single, _1, lambda, nu, log_z));
   }
};

struct Pcmp : public Worker {   
   // source vector and parameters
   const RVector<double> input;
   const double lambda;
   const double nu;
   const double log_z;
   
   // destination vector
   RVector<double> output;
   
   // initialize with source and destination
   Pcmp(const NumericVector input, const double lambda, const double nu, const double log_z, NumericVector output) 
      : input(input), lambda(lambda), nu(nu), log_z(log_z), output(output) {}
   
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::boost::bind(&pcmp_single, _1, lambda, nu, log_z));
   }
};

struct Qcmp : public Worker {   
   // source vector and parameters
   const RVector<double> input;
   const double lambda;
   const double nu;
   const double log_z;
   
   // destination vector
   RVector<double> output;
   
   // initialize with source and destination
   Qcmp(const NumericVector input, const double lambda, const double nu, const double log_z, NumericVector output) 
      : input(input), lambda(lambda), nu(nu), log_z(log_z), output(output) {}
   
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::boost::bind(&qcmp_single, _1, lambda, nu, log_z));
   }
};

#endif
