// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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

// [[Rcpp::export]]
NumericMatrix generate_attributes_cpp(int n, int m, double lambda) {
  // Use RcppArmadillo for efficient matrix operations
  arma::mat x;
  bool valid_matrix = false;

  while (!valid_matrix) {
    // Generate matrix of normal random values
    x = arma::randn(n, m) * std::sqrt(1.0/lambda);

    if (m == 1) {
      valid_matrix = true;
    } else {
      // Check if any option dominates another
      valid_matrix = true;
      arma::rowvec max_vals = arma::max(x, 0); // max for each column

      for (int i = 0; i < n; i++) {
        int count = 0;
        for (int j = 0; j < m; j++) {
          if (x(i,j) == max_vals(j)) count++;
        }
        if (count == m) {
          valid_matrix = false;
          break;
        }
      }
    }
  }

  // Convert to R matrix
  return Rcpp::wrap(x);
}

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

    // Update tracking
    fix_sequence(t) = current_fix;
    t++;
    thresh += delta;

    // Check termination conditions
    arma::vec opt_mean = mu * w;
    arma::vec opt_var = (1.0/prec) * w2;
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

    if(should_terminate || t >= max_steps) break;
  }

  // Prepare return values
  arma::uvec final_fix_sequence(t);
  for(int i = 0; i < t; i++) {
    final_fix_sequence(i) = fix_sequence(i);
  }

  // Calculate option values
  arma::vec opt_values = trial_x * w;
  int best_opt = opt_values.index_max();

  // Calculate fixation proportions
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
