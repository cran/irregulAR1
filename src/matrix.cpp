#include <RcppArmadillo.h>
#include <cmath>
#include "matrix.h"
// [[Rcpp::depends(RcppArmadillo)]]


// Covariance matrices ---------------------------------------------------------

arma::mat ar1_cov_consecutive(const int n,
                              const double rho,
                              const double sigma) {
  arma::mat A(n, n);
  double t;
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      t = static_cast<double>(i - j);
      A(i, j) = std::pow(rho, t) * pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

arma::mat ar1_cov_irregular(const arma::uvec& times,
                            const double rho,
                            const double sigma) {
  arma::mat A(times.size(), times.size());
  double t1, t2;
  for (int j = 0; j < times.size(); ++j) {
    t1 = static_cast<double>(times(j));
    for (int i = j + 1; i < times.size(); ++i) {
      t2 = static_cast<double>(times(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
                pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

arma::mat ar1_cross_cov(const arma::uvec& times1,
                        const arma::uvec& times2,
                        const double rho,
                        const double sigma) {
  int n = times1.size();
  int m = times2.size();
  arma::mat A(m, n);

  double t1, t2;
  for (int j = 0; j < n; ++j) {
    t1 = static_cast<double>(times1(j));
    for (int i = 0; i < m; ++i) {
      t2 = static_cast<double>(times2(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
        pow(sigma, 2) / (1 - std::pow(rho, 2));
    }
  }
  return A;
}

arma::mat ar1_cov_chol_irregular(const arma::uvec& times,
                                 const double rho,
                                 const double sigma) {
  return arma::chol(ar1_cov_irregular(times, rho, sigma));
}

// Precision matrices ----------------------------------------------------------

arma::sp_mat ar1_prec_consecutive(const int n,
                                  const double rho,
                                  const double sigma) {
  arma::sp_mat Q(n, n);
  Q(0, 0) = 1.0 / std::pow(sigma, 2.0);
  for (int i = 1; i < n; ++i) {
    Q(i, i - 1) = -rho / std::pow(sigma, 2.0);
    Q(i - 1, i) = -rho / std::pow(sigma, 2.0);
    Q(i, i)     = (1.0 + std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  }
  Q(n-1, n-1) = 1.0 / std::pow(sigma, 2.0);
  return Q;
}

arma::sp_mat ar1_prec_irregular(const arma::uvec& times,
                                const double rho,
                                const double sigma) {
  int n = times.size();
  arma::sp_mat Q(n, n);
  double a = (1.0 - std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  Q(0, 0) = a / (1.0 - std::pow(rho, 2.0 * (times(1) - times(0))));
  Q(n-1, n-1) = a / (1.0 - std::pow(rho, 2.0 * (times(n-1) - times(n-2))));
  double num, den;
  // Fill diagonal
  for (int i = 1; i < n - 1; ++i) {
    num = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i-1)));
    den = (1.0 - std::pow(rho, 2.0 * (times(i) - times(i-1))))
      * (1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i))));
    Q(i, i) = a * num / den;
  }
  // Fill first diagonals above and below main diagonal
  for (int i = 0; i < n - 1; ++i) {
    num = std::pow(rho, times(i+1) - times(i));
    den = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i)));
    Q(i + 1, i) = -a * num / den;
    Q(i, i + 1) = Q(i + 1, i);
  }
  return Q;
}

arma::sp_mat chol_tridiag_upper(const arma::sp_mat& Q) {
  int n = Q.n_rows;
  arma::sp_mat U(n, n);
  int p = 1;
  arma::vec v(n);
  for (int j = 0; j < n; ++j) {
    int lambda = std::min(j + p, n - 1);
    for (auto k = j, i = 0; k <= lambda; ++k, ++i) {
      v(k) = Q(k, j);
    }
    for (int k = std::max(0, j - p); k < j; ++k) {
      int i = std::min(k + p, n - 1);
      for (int h = 0; h <= i; ++h) {
        v(h) -= U(k, h) * U(k, j);
      }
    }
    U(j, arma::span(j, lambda)) = v.subvec(j, lambda).t() / std::sqrt(v(j));
  }
  return U;
}

arma::sp_mat ar1_prec_chol_irregular(const arma::uvec& times,
                                     const double rho,
                                     const double sigma) {
  return chol_tridiag_upper(ar1_prec_irregular(times, rho, sigma));
}

arma::vec band1_backsolve_vec(const arma::sp_mat& U,
                              const arma::vec& z) {
  int m = z.size();
  arma::vec v(m);
  v(m-1) = z(m-1) / U(m-1, m-1);
  for (int i = m - 2; i >= 0; --i) {
    v(i) = (z(i) - U(i, i + 1) * v(i+1)) / U(i, i);
  }
  return v;
}

arma::sp_mat band1_backsolve_mat(const arma::sp_mat& L,
                                 const arma::sp_mat& Q) {
  int m = L.n_cols;
  arma::sp_mat Y(m, m);

  // Fill diagonal and superdiagonal of first column
  Y(0, 0) = Q(0, 0) / L(0, 0);
  Y(1, 0) = (Q(1, 0) - L(1, 0) * Y(0, 0));

  // Fill sub-, main- and superdiagonals of all but first column
  for (int j = 1; j < m; ++j) {
    Y(j-1, j) =  Q(j-1, j) / L(j-1, j-1); // superdiag
    Y(j, j-1) = (Q(j, j-1) - L(j, j-1) * Y(j-1, j-1)) / L(j, j); // subdiag
    Y(j, j)   = (Q(j, j) - L(j, j-1) * Y(j-1, j)) / L(j, j); // diag
  }

  // Fill below subdiagonal
  for (int j = 0; j < m; ++j) {
    for (int i = j + 2; i < m; ++i) {
      Y(i, j) = -Y(i-1, j) * L(i, i-1) / L(i, i);
    }
  }
  return Y;
}

// Helper functions for derivatives of precision matrix Q

double dQ_corner(const double r, const double s,
                 const arma::uword t1, const arma::uword t2) {
  auto dt = t2 - t1;
  double a = std::pow(s, -2.0);
  double b = 1.0 - std::pow(r, 2.0 * dt);
  double term1 = dt * (1 - std::pow(r, 2.0))
               * std::pow(r, 2 * dt - 1.0) * std::pow(b, -2.0);
  double term2 = r / (1 - std::pow(r, 2 * dt));
  return 2 * a * (term1 - term2);
}

double dQ_diag(const double r, const double s,
               const arma::uword t1,
               const arma::uword t2,
               const arma::uword t3) {
  double a = std::pow(s, -2.0);
  auto dt1 = t2 - t1;
  auto dt2 = t3 - t2;
  auto dt3 = t3 - t1;
  auto den1 = 1.0 - std::pow(r, 2.0 * dt1);
  auto den2 = 1.0 - std::pow(r, 2.0 * dt2);
  double term1 = -r * (1 - std::pow(r, 2.0 * dt3)) / (den1 * den2);
  double term2 = (1 - std::pow(r, 2.0)) * dt1
               * std::pow(r, 2.0 * dt1 - 1.0) * (1 - std::pow(r, 2.0 * dt3))
               / (std::pow(den1, 2.0) * den2);
  double term3 = -(1.0 - std::pow(r, 2.0)) * dt3 * std::pow(r, 2.0 * dt3 - 1.0)
               / (den1 * den2);
  double term4 = (1.0 - std::pow(r, 2.0)) * dt2
               * (1.0 - std::pow(r, 2.0 * dt3)) * std::pow(r, 2.0 * dt2 - 1.0)
               / (den1 * std::pow(den2, 2.0));
  return 2 * a * (term1 + term2 + term3 + term4);
}

double dQ_offdiag(const double r, const double s,
                  const arma::uword t1,
                  const arma::uword t2) {
  double a = std::pow(s, -2.0);
  auto dt = t2 - t1;
  auto den = 1.0 - std::pow(r, 2.0 * dt);
  double term1 = 2.0 * std::pow(r, dt + 1.0) / den;
  double term2 = -(1.0 - std::pow(r, 2.0)) * dt * std::pow(r, dt - 1.0) / den;
  double term3 = -2.0 * (1.0 - std::pow(r, 2.0)) * dt
               * std::pow(r, 3 * dt - 1.0) * std::pow(den, -2.0);
  return a * (term1 + term2 + term3);
}

arma::sp_mat dprec_drho(const arma::uvec& times,
                        const double rho,
                        const double sigma) {
  int n = times.size();
  arma::sp_mat dQ(n, n);
  dQ(0, 0)     = dQ_corner(rho, sigma, times(0), times(1));
  dQ(n-1, n-1) = dQ_corner(rho, sigma, times(n-2), times(n-1));
  double den;
  // Fill diagonal
  for (int i = 1; i < n - 1; ++i) {
    dQ(i, i) = dQ_diag(rho, sigma, times(i-1), times(i), times(i+1));
  }
  // Fill first diagonals above and below main diagonal
  for (int i = 0; i < n - 1; ++i) {
    dQ(i + 1, i) = dQ_offdiag(rho, sigma, times(i), times(i+1));
    dQ(i, i + 1) = dQ(i + 1, i);
  }
  return dQ;
}

arma::sp_mat mult_U_band1U(const arma::sp_mat& A, const arma::sp_mat U) {
  auto m = U.n_cols;
  arma::sp_mat X(m, m);
  X(0, 0) = A(0, 0) * U(0, 0);
  for (int j = 1; j < m; ++j) {
    X(j, j) = A(j, j) * U(j, j);
    X(j-1, j) = A(j-1, j-1) * U(j-1, j) + A(j-1, j) * U(j, j);
  }
  return X;
}

arma::sp_mat dprecchol_drho(const arma::sp_mat& U, const arma::sp_mat& dQ) {
  arma::sp_mat B = band1_backsolve_mat(U.t(), dQ).t();
  arma::sp_mat A = band1_backsolve_mat(U.t(), B).t();
  for (int i = 0; i < A.n_cols; ++i) A(i, i) *= 0.5;
  return mult_U_band1U(A, U);
}
