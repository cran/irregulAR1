#ifndef MATRIX_H
#define MATRIX_H

#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

//' Covariance matrix for a stationary Gaussian AR(1) process, observed at
//' consecutive timepoints.
//'
//' Creates the covariance matrix of an AR(1) process with parameters \code{rho}
//' and \code{sigma}, observed at \code{n} consecutive time points. The process
//' is assumed to be in stationarity and to have Gaussian errors.
//' @param n An integer greater than or equal to 1.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{n} rows and \code{n} columns.
//' @export
//' @examples
//' n <- 5
//' rho <- 0.5
//' sigma <- 1
//' ar1_cov_consecutive(n, rho, sigma)
// [[Rcpp::export]]
arma::mat ar1_cov_consecutive(const int n,
                              const double rho,
                              const double sigma);

//' Covariance matrix for a stationary Gaussian AR(1) process, observed at
//' irregularly spaced time points.
//'
//' Creates the covariance matrix of an AR(1) process with parameters \code{rho}
//' and \code{sigma}, observed at the time points in the vector \code{times}.
//' The process is assumed to be in stationarity and to have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
//' @examples
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' ar1_cov_irregular(times, rho, sigma)
// [[Rcpp::export]]
arma::mat ar1_cov_irregular(const arma::uvec& times,
                            const double rho,
                            const double sigma);

//' Cross-covariance matrix of a stationary Gaussian AR(1) process.
//'
//' Creates the cross-covariance matrix of an AR(1) process with parameters
//' \code{rho} and \code{sigma}, observed at (positive) integer times
//' \code{times1} and \code{times2}, which may be irregularly spaced. The
//' process is assumed to be in stationarity and to have Gaussian errors.
//' @param times1 An vector of positive integers, preferably ordered.
//' @param times2 An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{length(times2)} rows and \code{length(times1)}
//'   columns.
//' @export
//' @examples
//' times1 <- c(1, 3, 6)
//' times2 <- c(2, 4, 8:9)
//' rho <- 0.5
//' sigma <- 1
//' ar1_cross_cov(times1, times2, rho, sigma)
// [[Rcpp::export]]
arma::mat ar1_cross_cov(const arma::uvec& times1,
                        const arma::uvec& times2,
                        const double rho,
                        const double sigma);

//' Upper triangular Cholesky decomposition for a stationary Gaussian AR(1)
//' process covariance matrix, observed at irregularly spaced time points.
//'
//' Creates the upper Cholesky triangle of the covariance matrix of an AR(1)
//' process with parameters \code{rho} and \code{sigma}, observed at the time
//' points in the vector \code{times}. The process is assumed to be in
//' stationarity and to have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
//' @examples
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' ar1_cov_chol_irregular(times, rho, sigma)
// [[Rcpp::export]]
arma::mat ar1_cov_chol_irregular(const arma::uvec& times,
                                 const double rho,
                                 const double sigma);

//' Sparse precision matrix for a stationary Gaussian AR(1) process, observed at
//' consecutive timepoints.
//'
//' Creates the precision (inverse covariance) matrix of an AR(1) process with
//' parameters \code{rho} and \code{sigma}, observed at \code{n} consecutive
//' time points. The process is assumed to be in stationarity and to have
//' Gaussian errors. The matrix is a tridiagonal band matrix and thus sparse.
//' @param n An integer greater than or equal to 1.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{n} rows and \code{n} columns.
//' @export
//' @examples
//' library(Matrix)
//' n <- 5
//' rho <- 0.5
//' sigma <- 1
//' ar1_prec_consecutive(n, rho, sigma)
// [[Rcpp::export]]
arma::sp_mat ar1_prec_consecutive(const int n,
                                  const double rho,
                                  const double sigma);

//' Precision matrix for a stationary Gaussian AR(1) process, observed at
//' irregularly spaced time points.
//'
//' Creates the precision (inverse covariance) matrix of an AR(1) process with
//' parameters \code{rho} and \code{sigma}, observed at the time points in the
//' vector \code{times}. The process is assumed to be in stationarity and to
//' have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
//' @examples
//' library(Matrix)
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' ar1_prec_irregular(times, rho, sigma)
// [[Rcpp::export]]
arma::sp_mat ar1_prec_irregular(const arma::uvec& times,
                                const double rho,
                                const double sigma);

//' Upper Cholesky decomposition of a tridiagonal matrix.
//'
//' Creates the lower Cholesky decomposition of a tridiagonal matrix. The
//' decomposition will be a sparse lower triangular matrix with non-zero
//' elements only on the main diagonal and the diagonal below it.
//' @param Q A square tridiagonal matrix.
//' @return A sparse square matrix with the same size as the input matrix.
//' @export
//' @examples
//' library(Matrix)
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' Q <- ar1_prec_irregular(times, rho, sigma)
//' chol_tridiag_upper(Q)
// [[Rcpp::export]]
arma::sp_mat chol_tridiag_upper(const arma::sp_mat& Q);

//' Upper Cholesky triangle of the precision matrix of a stationary Gaussian
//' AR(1) process, observed at irregularly spaced time points.
//'
//' Creates the upper triangular Cholesky decomposition of the precision matrix
//' of an AR(1) process with parameters \code{rho} and \code{sigma}, observed at
//' the time points in the vector \code{times}. The process is assumed to be in
//' stationarity and to have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A sparse square matrix with \code{length(times)} rows.
//' @export
//' @examples
//' library(Matrix)
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' ar1_prec_chol_irregular(times, rho, sigma)
// [[Rcpp::export]]
arma::sp_mat ar1_prec_chol_irregular(const arma::uvec& times,
                                     const double rho,
                                     const double sigma);


//' Backsolve with band 1 upper Cholesky.
//'
//' Backsolve with band 1 upper Cholesky.
//' @param U An upper triangular square matrix with non-zero entries only on the
//'   main diagonal and the first superdiagonal.
//' @param z A vector with as many elements as the number of rows of U.
//' @return A vector.
//' @keywords internal
// [[Rcpp::export]]
arma::vec band1_backsolve_vec(const arma::sp_mat& U,
                              const arma::vec& z);

//' Backward substitution with band 1 lower Cholesky triangle and tridiagonal
//' RHS.
//'
//' Backward substitution with band 1 lower Cholesky triangle and tridiagonal
//' matrix on the right hand side.
//' @param L A lower triangular square matrix with non-zero entries only on the
//'   main diagonal and the first subdiagonal.
//' @param Q A tridiagonal matrix with the same dimensions as L.
//' @return A matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat band1_backsolve_mat(const arma::sp_mat& L,
                                 const arma::sp_mat& Q);

//' Derivative of the precision matrix for a stationary Gaussian AR(1) process.
//'
//' Creates the derivate of the precision matrix of an AR(1) process with
//' respect to the parameter \code{rho}. The process has been observed at the
//' time points in the vector \code{times} and is assumed to be in stationarity,
//' and to have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A sparse square matrix with \code{length(times)} rows.
//' @export
//' @examples
//' library(Matrix)
//' times <- c(1, 4:5, 7)
//' rho <- 0.5
//' sigma <- 1
//' dprec_drho(times, rho, sigma)
// [[Rcpp::export]]
arma::sp_mat dprec_drho(const arma::uvec& times,
                        const double rho,
                        const double sigma);

//' Multiply an upper triangular matrix with a band 1 upper triangular matrix.
//'
//' Multiply an upper triangular matrix with a band 1 upper triangular matrix.
//' @param A A sparse upper triangular matrix.
//' @param U A sparse band 1 upper triangular matrix of the same dimensions as
//'   \code{A}.
//' @return A sparse band 1 upper triangular matrix.
// [[Rcpp::export]]
arma::sp_mat mult_U_band1U(const arma::sp_mat& A, const arma::sp_mat U);

//' Derivative of the upper Cholesky triangle of the precision matrix of a
//' stationary Gaussian AR(1) process.
//'
//' Creates the derivate of the upper Cholesky triangle of the precision matrix
//' of an AR(1) process with respect to the parameter \code{rho}.
//' @param U The upper Cholesky triangle of the precision matrix \code{Q} of the
//'   AR(1) process.
//' @param dQ The derivative of the precision matrix \code{Q} with respect to
//'   the correlation parameter \code{rho}.
//' @return A band 1 upper triangular matrix of the same dimensions as \code{U}.
//' @export
//' @examples
//' library(Matrix)
//' t <- c(1, 3:4, 6, 8)
//' r <- 0.5
//' s <- 1
//' U <- ar1_prec_chol_irregular(t, r, s)
//' dQ <- dprec_drho(t, r, s)
//' (dU <- dprecchol_drho(U, dQ))
// [[Rcpp::export]]
arma::sp_mat dprecchol_drho(const arma::sp_mat& U, const arma::sp_mat& dQ);

#endif
