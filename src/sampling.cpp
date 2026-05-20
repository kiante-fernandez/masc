#include "masc_rcpp.h"

// In-place Kalman update of one option's multivariate-normal belief, in
// covariance form: innovation variance S = Sigma(j,j) + sigma_s^2, gain
// K = Sigma[,j] / S, mean shift mu += K * (sample - mu(j)), covariance shrink
// Sigma -= S * K K'. Observing attribute j therefore moves the belief about
// every attribute correlated with j. Single source of truth shared by the
// sampling loop and the exported masc_kalman_update_cpp() used in tests; `j` is
// 0-based. When Sigma is diagonal K is zero off attribute j, so this reduces
// exactly to the univariate normal-normal precision-weighted update.
static void masc_kalman_update_inplace(arma::vec& mu, arma::mat& Sigma, int j,
                                       double sample, double sigma_s_sq) {
  double S = Sigma(j, j) + sigma_s_sq;       // innovation variance
  arma::vec K = Sigma.col(j) / S;            // Kalman gain (vector)
  double innovation = sample - mu(j);
  mu += K * innovation;                       // belief spread to all attributes
  Sigma = Sigma - S * (K * K.t());           // covariance shrink
  Sigma = 0.5 * (Sigma + Sigma.t());         // re-symmetrize
}

// [[Rcpp::export]]
List masc_kalman_update_cpp(arma::vec mu, arma::mat Sigma, int j,
                            double sample, double sigma_s_sq) {
  // `j` is 1-based when called from R
  masc_kalman_update_inplace(mu, Sigma, j - 1, sample, sigma_s_sq);
  return List::create(Named("mu") = mu, Named("Sigma") = Sigma);
}

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
                       int n_attributes,
                       Rcpp::Nullable<Rcpp::NumericMatrix> Sigma_belief = R_NilValue) {

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

  // ----------------------------------------------------------------------
  // Multivariate (MASC-C) belief update branch.
  //
  // Active only when a non-diagonal belief covariance is supplied. Beliefs for
  // each option are a multivariate normal tracked by a mean row in `mu` and a
  // covariance matrix in `Sigma_list`, updated via the Kalman filter so that
  // observing one attribute spreads information to correlated attributes.
  //
  // When Sigma_belief is NULL or diagonal we fall through to the original
  // univariate code below, which keeps that path byte-for-byte identical
  // (a diagonal covariance makes the multivariate update reduce to it exactly).
  // The myopic search rule is reused unchanged: it receives the diagonal
  // precisions 1/Sigma_jj extracted from each option's current covariance.
  // ----------------------------------------------------------------------
  bool multivariate = false;
  arma::mat Sigma_b;
  if (Sigma_belief.isNotNull()) {
    Rcpp::NumericMatrix Sb(Sigma_belief);
    if (Sb.nrow() != n_attributes || Sb.ncol() != n_attributes) {
      stop("Error: Sigma_belief must be %d x %d", n_attributes, n_attributes);
    }
    Sigma_b = arma::mat(Sb.begin(), n_attributes, n_attributes, true);
    double off_diag_max = 0.0;
    for (int a = 0; a < n_attributes; ++a) {
      for (int b = 0; b < n_attributes; ++b) {
        if (a != b) off_diag_max = std::max(off_diag_max, std::abs(Sigma_b(a, b)));
      }
    }
    multivariate = (off_diag_max > 1e-10);
  }

  if (multivariate) {
    // Initialize belief distributions: mu_i = 0, Sigma_i = Sigma_belief / lambda
    arma::mat mu(n_options, n_attributes, arma::fill::zeros);
    const arma::mat Sigma0 = Sigma_b / lambda;
    std::vector<arma::mat> Sigma_list(n_options, Sigma0);

    // Marginal precisions 1/Sigma_jj feed the (unchanged) myopic search rule.
    // Built once, then refreshed only for the option whose belief just changed.
    arma::mat prec(n_options, n_attributes);
    for (int i = 0; i < n_options; ++i)
      for (int j = 0; j < n_attributes; ++j)
        prec(i, j) = 1.0 / Sigma_list[i](j, j);

    // Loop-invariant search-rule inputs, wrapped once.
    const NumericVector w_r = wrap(w), w2_r = wrap(w2), sp_r = wrap(sp);

    int t = 0;
    double thresh = theta;
    arma::uvec fix_sequence(max_steps, arma::fill::zeros);

    while (t < max_steps) {
      // Same search rule as the univariate model
      NumericVector trans_probs = MASC_SearchRule_myopic_cpp(
        n_options, n_attributes, w_r, w2_r, sp_r,
        thresh, alpha, wrap(prec), wrap(mu)
      );
      arma::vec trans_vec(trans_probs.begin(), trans_probs.size(), false);

      // Sample current fixation
      double u = R::runif(0, 1);
      double cumsum = 0;
      int current_fix = 0;
      for (int i = 0; i < (int) trans_vec.n_elem; i++) {
        cumsum += trans_vec(i);
        if (u <= cumsum) {
          current_fix = i;
          break;
        }
      }

      int j_fix = current_fix / n_options;  // Attribute index (column)
      int i_fix = current_fix % n_options;  // Option index (row)

      // --- Sample ---
      double sample = trial_x(i_fix, j_fix) + R::rnorm(0, sigma(j_fix));

      // --- Kalman update of the fixated option's belief (covariance form) ---
      arma::vec mu_i = mu.row(i_fix).t();
      masc_kalman_update_inplace(mu_i, Sigma_list[i_fix], j_fix, sample,
                                 sigma(j_fix) * sigma(j_fix));
      mu.row(i_fix) = mu_i.t();
      for (int j = 0; j < n_attributes; ++j)            // refresh changed option only
        prec(i_fix, j) = 1.0 / Sigma_list[i_fix](j, j);

      // Termination: same decision rule as the univariate path, but option
      // variance is the full quadratic form w' Sigma_i w (includes off-diagonal
      // covariance), not the diagonal sum.
      arma::vec opt_mean = mu * w;
      arma::vec opt_var = zeros(n_options);
      for(int i = 0; i < n_options; i++) {
        opt_var(i) = arma::as_scalar(w.t() * Sigma_list[i] * w);
      }

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

    arma::vec opt_values = trial_x * w;
    int best_opt = opt_values.index_max();

    arma::vec prop_fix_opt(n_options, arma::fill::zeros);
    arma::vec prop_fix_att(n_attributes, arma::fill::zeros);
    if(t > 0) {
      for(int i = 0; i < t; ++i) {
        int fix = fix_sequence[i];
        prop_fix_opt(fix % n_options) += 1.0/t;
        prop_fix_att(fix / n_options) += 1.0/t;
      }
    }

    return List::create(
      Named("response") = mu * w,
      Named("best_option") = best_opt + 1,
      Named("rt") = t,
      Named("fix_sequence") = final_fix_sequence + 1,
      Named("prop_fix_opt") = prop_fix_opt,
      Named("prop_fix_att") = prop_fix_att
    );
  }

  // ===================== Univariate (original MASC) path =====================
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
