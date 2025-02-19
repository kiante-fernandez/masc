#include "masc_rcpp.h"

// [[Rcpp::export]]
List rMASC_sampling_cpp(const arma::mat& trial_x,
                       const arma::vec& w,
                       double sigma,
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
  // Pre-compute squared weights and sampling precision
  arma::vec w2 = w % w;
  arma::vec sp = arma::vec(n_attributes).fill(1.0/(sigma*sigma));

  // Initialize belief distributions
  arma::mat prec = arma::mat(n_options, n_attributes).fill(lambda);
  arma::mat mu = arma::mat(n_options, n_attributes).zeros();

  // Initialize tracking variables
  int t = 0;
  double thresh = theta;
  arma::uvec fix_sequence(max_steps);

  // Main decision loop
  while(true) {
    // Get transition probabilities
    NumericVector trans_probs = MASC_SearchRule_myopic_cpp(
      n_options, n_attributes,
      wrap(w), wrap(w2), wrap(sp),
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

    // Calculate attribute being fixated
    int j_fix = current_fix / n_options;

    // Sample and update beliefs
    double current_sample = trial_x(current_fix) + R::rnorm(0, sigma);
    double new_prec_val = prec(current_fix) + sp(j_fix);
    mu(current_fix) = (current_sample * sp(j_fix) +
      mu(current_fix) * prec(current_fix)) / new_prec_val;
    prec(current_fix) = new_prec_val;

    // Check termination conditions
    arma::vec opt_mean = mu * w;
    arma::vec opt_var = zeros(n_options);
    for(int i = 0; i < n_options; i++) {
      for(int j = 0; j < n_attributes; j++) {
        opt_var(i) += (1.0/prec(i,j)) * w2(j);
      }
    }

    int current_best = opt_mean.index_max();

    bool should_terminate = true;
    for(int i = 0; i < n_options; i++) {
      if(i != current_best) {
        double z_score = (opt_mean(current_best) - opt_mean(i)) /
          std::sqrt(opt_var(current_best) + opt_var(i));
        if(R::pnorm(0, z_score, 1, 1, 0) >= thresh) {
          should_terminate = false;
          break;
        }
      }
    }
        // Update tracking
    fix_sequence(t) = current_fix;
    t++;
    thresh += delta;
    if(should_terminate || t >= max_steps) break;
  }

  // Prepare return values
  arma::uvec final_fix_sequence(t);
  for(int i = 0; i < t; i++) {
    final_fix_sequence(i) = fix_sequence(i);
  }

  // Calculate option values and fixation proportions
  arma::vec opt_values = trial_x * w;
  int best_opt = opt_values.index_max();

  arma::vec prop_fix_opt(n_options, arma::fill::zeros);
  arma::vec prop_fix_att(n_attributes, arma::fill::zeros);

  for(int i = 0; i < t; i++) {
    int opt = fix_sequence(i) % n_options;
    int att = fix_sequence(i) / n_options;
    prop_fix_opt(opt) += 1.0/t;
    prop_fix_att(att) += 1.0/t;
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
