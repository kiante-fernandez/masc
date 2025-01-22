#include "masc_rcpp.h"

// [[Rcpp::export]]
NumericVector MASC_SearchRule_myopic_cpp(
    int n, int m,
    NumericVector w,
    NumericVector w2,
    NumericVector sp,
    double thresh,
    double alpha,
    NumericMatrix prec,
    NumericMatrix mu) {

  // Convert R objects to Armadillo
  arma::vec w_vec(w.begin(), w.length(), false);
  arma::vec w2_vec(w2.begin(), w2.length(), false);
  arma::vec sp_vec(sp.begin(), sp.length(), false);
  arma::mat prec_mat(prec.begin(), n, m, false);
  arma::mat mu_mat(mu.begin(), n, m, false);

  // Precompute shared terms
  arma::mat new_prec = prec_mat.each_row() + sp_vec.t();

  arma::vec opt_means = mu_mat * w_vec;
  arma::vec opt_vars_old(n, arma::fill::zeros);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      opt_vars_old(i) += w2_vec(j) / prec_mat(i,j);
    }
  }

  arma::mat myopic_score(n, m, arma::fill::zeros);

  for(int i = 0; i < n; ++i) {
    arma::uvec not_i = arma::find(arma::linspace<arma::uvec>(0, n-1, n) != i);

    for(int j = 0; j < m; ++j) {
      // Calculate new variance
      double opt_var_new = w2_vec(j) / new_prec(i,j);
      if(m > 1) {
        arma::uvec not_j = arma::find(arma::linspace<arma::uvec>(0, m-1, m) != j);
        for(arma::uword k = 0; k < not_j.n_elem; ++k) {
          opt_var_new += w2_vec(not_j(k)) / prec_mat(i, not_j(k));
        }
      }

      // MATLAB-style variance calculation
      arma::vec var_term = opt_var_new + opt_vars_old.elem(not_i);

      // Threshold computation
      arma::vec opt_mean_thresh(not_i.n_elem);
      for(arma::uword k = 0; k < not_i.n_elem; ++k) {
        double sd = std::sqrt(var_term(k));
        opt_mean_thresh(k) = opt_means(not_i(k)) - R::qnorm(thresh, 0.0, sd, 1, 0);
      }

      // Calculate other_terms ONCE PER ATTRIBUTE (FIXED)
      double other_terms = 0.0;
      if(m > 1) {
        arma::uvec not_j = arma::find(arma::linspace<arma::uvec>(0, m-1, m) != j);
        for(arma::uword l = 0; l < not_j.n_elem; ++l) {
          other_terms += mu_mat(i, not_j(l)) * w_vec(not_j(l));
        }
      }

      // MATLAB-style threshold calculation
      arma::vec threshold_terms(not_i.n_elem);
      for(arma::uword k = 0; k < not_i.n_elem; ++k) {
        double term = (new_prec(i,j)/w_vec(j)) * (opt_mean_thresh(k) - other_terms);
        threshold_terms(k) = (term - prec_mat(i,j)*mu_mat(i,j)) / sp_vec(j);
      }

      double sample_thresh = arma::max(threshold_terms);

      // MATLAB-style CDF calculation
      double att_sd = std::sqrt(1.0 / prec_mat(i,j));
      myopic_score(i,j) = R::pnorm(mu_mat(i,j), sample_thresh, att_sd, 1, 0);
    }
  }

  // MATLAB-style normalization
  myopic_score.replace(0.0, arma::datum::eps);
  myopic_score /= arma::accu(myopic_score);

  arma::mat transition_prob = arma::exp(alpha * myopic_score);
  transition_prob /= arma::accu(transition_prob);

  return NumericVector(wrap(arma::vectorise(transition_prob)));
}
