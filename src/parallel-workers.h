#ifndef PARWORK
#define PARWORK

#include <Rcpp.h>
#include <math.h>
#include <RcppParallel.h>
#include <boost/bind.hpp>
#include "compoisson.h"

using namespace RcppParallel;


struct Dcom : public Worker {   
   // source vector and parameters
   const RVector<double> input;
   const double lambda;
   const double nu;
   const double log_z;
   
   // destination vector
   RVector<double> output;
   
   // initialize with source and destination
   Dcom(const NumericVector input, const double lambda, const double nu, const double log_z, NumericVector output) 
      : input(input), lambda(lambda), nu(nu), log_z(log_z), output(output) {}
   
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::boost::bind(&dcom_single, _1, lambda, nu, log_z));
   }
};

#endif
