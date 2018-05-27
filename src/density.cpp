#include <cmath>
#include "density.h"
#include "matrix.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double ar1_lpdf_cpp(const arma::vec& x,
                    const arma::vec& mu,
                    const arma::uvec& times,
                    const double rho,
                    const double sigma) {
  arma::sp_mat Q = ar1_prec_irregular(times, rho, sigma);
  arma::sp_mat U = chol_tridiag_upper(Q);

  double diag_sum = 0.0;
  for (int i = 0; i < U.n_rows; ++i) diag_sum += std::log(U(i, i));

  arma::vec z = U * (x - mu);
  double q = arma::dot(z, z);

  double c = -0.5 * U.n_rows * std::log(2.0 * arma::datum::pi);

  return c + diag_sum - q * 0.5;
}
