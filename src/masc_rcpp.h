#ifndef MASC_RCPP_H
#define MASC_RCPP_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Function declarations
NumericVector MASC_SearchRule_myopic_cpp(
    int n, int m,
    NumericVector w,
    NumericVector w2,
    NumericVector sp,
    double thresh,
    double alpha,
    NumericMatrix prec,
    NumericMatrix mu);

NumericMatrix generate_attributes_cpp(int n, int m, double lambda);

List rMASC_sampling_cpp(
    const arma::mat& trial_x,
    const arma::vec& w,
    double sigma,
    double alpha,
    double delta,
    double theta,
    double lambda,
    int max_steps,
    int n_options,
    int n_attributes);

#endif
