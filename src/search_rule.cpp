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
  vec w_vec(w.begin(), w.length(), false);
  vec w2_vec(w2.begin(), w2.length(), false);
  vec sp_vec(sp.begin(), sp.length(), false);
  mat prec_mat(prec.begin(), n, m, false);
  mat mu_mat(mu.begin(), n, m, false);

  // Precompute shared terms
  mat new_prec(n, m);
  for(int j = 0; j < m; j++) {
    new_prec.col(j) = prec_mat.col(j) + sp_vec(j);
  }

  vec opt_means = mu_mat * w_vec;
  vec opt_vars_old = zeros(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      opt_vars_old(i) += (1.0 / prec_mat(i,j)) * w2_vec(j);
    }
  }

  // Initialize score matrix
  mat myopic_score = zeros(n, m);

  // Calculate myopic score for each option-attribute pair
  for(int i = 0; i < n; i++) {
    uvec not_i = find(linspace<uvec>(0, n-1, n) != i);

    for(int j = 0; j < m; j++) {
      // Calculate new option variance
      double opt_var_new;
      if(m > 1) {
        uvec not_j = find(linspace<uvec>(0, m-1, m) != j);
        opt_var_new = (1.0/new_prec(i,j)) * w2_vec(j);
        for(uword k = 0; k < not_j.n_elem; k++) {
          opt_var_new += (1.0/prec_mat(i,not_j(k))) * w2_vec(not_j(k));
        }
      } else {
        opt_var_new = (1.0/new_prec(i,j)) * w2_vec(j);
      }

      // Safe variance calculation and threshold computation
      vec var_term = opt_var_new + opt_vars_old.elem(not_i);
      var_term = max(var_term, vec(var_term.n_elem).fill(datum::eps));

      vec opt_mean_thresh(not_i.n_elem);
      for(uword k = 0; k < not_i.n_elem; k++) {
        opt_mean_thresh(k) = opt_means(not_i(k)) -
          R::qnorm(thresh, 0.0, std::sqrt(var_term(k)), 1, 0);
      }

      // Calculate threshold
      double sample_thresh = 0.0;
      if(w_vec(j) > 0) {
        if(m > 1) {
          uvec not_j = find(linspace<uvec>(0, m-1, m) != j);
          double other_terms = 0.0;
          for(uword k = 0; k < not_j.n_elem; k++) {
            other_terms += mu_mat(i,not_j(k)) * w_vec(not_j(k));
          }
          sample_thresh = (new_prec(i,j)/w_vec(j)) *
            (mean(opt_mean_thresh) - other_terms);
          sample_thresh = (sample_thresh - prec_mat(i,j)*mu_mat(i,j))/sp_vec(j);
        } else {
          sample_thresh = (new_prec(i,j)/w_vec(j) * mean(opt_mean_thresh) -
            prec_mat(i,j)*mu_mat(i,j))/sp_vec(j);
        }
      }

      // Update score
      double att_sd = std::sqrt(1.0/prec_mat(i,j));
      myopic_score(i,j) = R::pnorm(mu_mat(i,j), std::max(sample_thresh, 0.0), att_sd, 1, 0);
    }
  }

  // Handle zero scores and normalize
  myopic_score.elem(find(myopic_score == 0)).fill(datum::eps);
  myopic_score = myopic_score / accu(myopic_score);

  // Apply search sensitivity and return probabilities
  mat transition_prob = exp(alpha * myopic_score);
  transition_prob = transition_prob / accu(transition_prob);

  // Convert to R vector
  NumericVector result = wrap(vectorise(transition_prob));
  return result;
}
