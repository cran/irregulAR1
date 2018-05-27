#ifndef DENSITY_H
#define DENSITY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Evaluate the log-density of a stationary Gaussian AR(1) process.
//'
//' Evaluate the log-density of a stationary Gaussian AR(1) process, observed at
//' times \code{times} taking values \code{x}.
//' @param x A vector of observed values.
//' @param mu A vector of expected values.
//' @param times A vector of the time points of observation.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A scalar, the log density.
//' @keywords internal
// [[Rcpp::export]]
double ar1_lpdf_cpp(const arma::vec& x,
                    const arma::vec& mu,
                    const arma::uvec& times,
                    const double rho,
                    const double sigma);

#endif
