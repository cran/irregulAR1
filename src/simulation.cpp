#ifndef SIMULATION_H
#define SIMULATION_H

#include <RcppArmadillo.h>
#include <cmath>
#include "simulation.h"
#include "matrix.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec ar1_sim_cpp(const int n,
                      const double rho,
                      const double sigma) {
  arma::vec x(n);
  x[0] = R::rnorm(0, 1) * sigma / std::sqrt(1 - std::pow(rho, 2));
  for (auto i = 1; i < n; ++i) {
    x[i] = rho * x[i - 1] + sigma * R::rnorm(0, 1);
  }
  return x;
}

arma::vec ar1_sim_irregular_cpp(const arma::uvec& times,
                                const double rho,
                                const double sigma) {
  arma::sp_mat Q = ar1_prec_irregular(times, rho, sigma);
  arma::vec z = Rcpp::as<arma::vec>(Rcpp::rnorm(Q.n_cols, 0.0, 1.0));
  arma::sp_mat U = chol_tridiag_upper(Q);
  return band1_backsolve_vec(U, z);
}

arma::vec ar1_sim_conditional_cpp(const arma::uvec& pred_times,
                                  const arma::vec& mu_pred,
                                  const arma::vec& x_obs,
                                  const arma::uvec& obs_times,
                                  const arma::vec& mu_obs,
                                  const double rho,
                                  const double sigma) {
  // Use sort_index on (pred_times, obs_times)
  arma::uvec all_times = arma::join_cols(pred_times, obs_times);
  arma::uvec sorted_times = arma::sort_index(all_times);

  // Create Q for the sorted vector
  arma::sp_mat Q_o = ar1_prec_irregular(arma::sort(obs_times), rho, sigma);

  // Draw x_all with this Q; separate into x_pred_star and x_obs_star
  arma::vec x_all = ar1_sim_irregular_cpp(arma::sort(all_times), rho, sigma);
  x_all = x_all(sorted_times);
  x_all = x_all + arma::join_cols(mu_pred, mu_obs);

  // Create Sigma_po
  arma::mat Sigma_po = ar1_cross_cov(obs_times, pred_times, rho, sigma);


  // return x_pred_star + Sigma_po * Q * (x_obs - x_obs_star)
  return x_all.head(pred_times.n_elem)
       + Sigma_po * Q_o * (x_obs - x_all.tail(obs_times.n_elem));
}

#endif
