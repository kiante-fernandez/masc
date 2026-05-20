#include "masc_rcpp.h"

// [[Rcpp::export]]
NumericMatrix generate_attributes_cpp(int n, int m, double lambda,
                                      Rcpp::Nullable<Rcpp::NumericMatrix> Sigma = R_NilValue) {
  // Use RcppArmadillo for efficient matrix operations
  arma::mat x;
  bool valid_matrix = false;

  // Optional correlation structure. When Sigma is NULL the draw is identical to
  // the original independent-attribute generator (preserving seeded behaviour):
  //   x = randn(n, m) * sqrt(1/lambda)
  // When supplied, rows are drawn from N(0, Sigma/lambda) via a Cholesky factor:
  //   x = randn(n, m) * R, where R'R = Sigma/lambda (upper-triangular arma::chol).
  bool use_sigma = Sigma.isNotNull();
  arma::mat chol_factor;
  if (use_sigma) {
    Rcpp::NumericMatrix S(Sigma);
    if (S.nrow() != m || S.ncol() != m) {
      stop("Sigma must be an m x m matrix (m = number of attributes)");
    }
    arma::mat Sigma_scaled(S.begin(), m, m, true);
    Sigma_scaled /= lambda;
    chol_factor = arma::chol(Sigma_scaled);  // upper-triangular R with R'R = Sigma_scaled
  }

  while (!valid_matrix) {
    // Generate matrix of normal random values
    if (use_sigma) {
      x = arma::randn(n, m) * chol_factor;
    } else {
      x = arma::randn(n, m) * std::sqrt(1.0/lambda);
    }

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
