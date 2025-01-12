#include "masc_rcpp.h"

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
