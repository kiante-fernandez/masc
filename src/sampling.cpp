#include "masc_rcpp.h"

// [[Rcpp::export]]
List rMASC_sampling_cpp(const arma::mat& trial_x,
                       const arma::vec& w,
                       const arma::vec& sigma,  // Now a vector!
                       double alpha,
                       double delta,
                       double theta,
                       double lambda,
                       int max_steps,
                       int n_options,
                       int n_attributes) {

  // Check dimensions
  if (trial_x.n_rows != n_options || trial_x.n_cols != n_attributes) {
    stop("Error: trial_x dimensions (%d, %d) do not match expected dimensions (%d, %d)",
         trial_x.n_rows, trial_x.n_cols, n_options, n_attributes);
  }
  if (w.n_elem != n_attributes) {
    stop("Error: weights vector length (%d) does not match number of attributes (%d)",
         w.n_elem, n_attributes);
  }
  if (sigma.n_elem != n_attributes) {
    stop("Error: sigma vector length (%d) does not match number of attributes (%d)",
         sigma.n_elem, n_attributes);
  }
  // Pre-compute squared weights and sampling precision
  const arma::vec w2 = arma::square(w);
  const arma::vec sp = 1.0 / arma::square(sigma);
  const int n_pairs = n_options * n_attributes;

  // Initialize belief distributions
  arma::mat prec(n_options, n_attributes, arma::fill::value(lambda));
  arma::mat mu(n_options, n_attributes, arma::fill::zeros);

  // Initialize tracking variables
  int t = 0;
  double thresh = theta;
  arma::uvec fix_sequence(max_steps, arma::fill::zeros);

  // Main decision loop
  while(t < max_steps) {
    // Get transition probabilities
    NumericVector trans_probs = MASC_SearchRule_myopic_cpp(
      n_options, n_attributes, wrap(w), wrap(w2), wrap(sp),
      thresh, alpha, wrap(prec), wrap(mu)
    );

    // Convert to arma vector for sampling
    arma::vec trans_vec(trans_probs.begin(), trans_probs.size(), false);
    // Sample current fixation
    double u = R::runif(0, 1);
    double cumsum = 0;
    int current_fix = 0;
    for(int i = 0; i < trans_vec.n_elem; i++) {
      cumsum += trans_vec(i);
      if(u <= cumsum) {
        current_fix = i;
        break;
      }
    }

    // Calculate attribute (j_fix) and option (i_fix) indices being fixated
    // These are needed for accessing trial_x and sigma
    int j_fix = current_fix / n_options;  // Attribute index (column)
    int i_fix = current_fix % n_options;  // Option index (row)

    // --- Sample ---
    // Access trial_x using (row, col) indices
    double sample = trial_x(i_fix, j_fix) + R::rnorm(0, sigma(j_fix));
    double old_prec = prec(i_fix, j_fix);
    prec(i_fix, j_fix) += sp(j_fix);
    mu(i_fix, j_fix) = (mu(i_fix, j_fix) * old_prec + sample * sp(j_fix)) / prec(i_fix, j_fix);

    // Check termination conditions
    arma::vec opt_mean = mu * w;
    arma::vec opt_var = zeros(n_options);
    for(int i = 0; i < n_options; i++) {
      for(int j = 0; j < n_attributes; j++) {
        opt_var(i) += (1.0/prec(i,j)) * w2(j);
      }
    }

    // Check for unique maximum (like MATLAB's sum(currentBest)==1)
    arma::uvec max_indices = arma::find(opt_mean == arma::max(opt_mean));
    bool unique_max = (max_indices.n_elem == 1);

    bool should_terminate = false;
    if (unique_max) {
      int current_best = max_indices(0);
      should_terminate = true;

      for(int i = 0; i < n_options; i++) {
        if(i != current_best) {
          double mu_diff = opt_mean(current_best) - opt_mean(i);
          double var_diff = opt_var(current_best) + opt_var(i);
          double sigma_diff = std::sqrt(var_diff);

          double prob = R::pnorm(0, mu_diff, sigma_diff, 1, 0);
          if(prob > thresh) {
            should_terminate = false;
            break;
          }
        }
      }
    }
        // Update tracking
    fix_sequence(t) = current_fix;
    t++;

    thresh += delta;

    if(should_terminate) break;
  }

  // Prepare return values
  arma::uvec final_fix_sequence(t);

  for(int i = 0; i < t; i++) {
    final_fix_sequence(i) = fix_sequence(i);
  }

  // Calculate option values and fixation proportions
  arma::vec opt_values = trial_x * w;
  int best_opt = opt_values.index_max();

  // Calculate proportions
  arma::vec prop_fix_opt(n_options, arma::fill::zeros);
  arma::vec prop_fix_att(n_attributes, arma::fill::zeros);
  if(t > 0) {
    for(int i = 0; i < t; ++i) {
      int fix = fix_sequence[i];
      prop_fix_opt(fix % n_options) += 1.0/t;
      prop_fix_att(fix / n_options) += 1.0/t;
    }
  }

  // Return results
  return List::create(
    Named("response") = mu * w,
    Named("best_option") = best_opt + 1,  // Convert to 1-based indexing
    Named("rt") = t,
    Named("fix_sequence") = final_fix_sequence + 1,  // Convert to 1-based
    Named("prop_fix_opt") = prop_fix_opt,
    Named("prop_fix_att") = prop_fix_att
  );
}
