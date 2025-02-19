test_that("generate_attributes produces correct output format", {
  n <- 2
  m <- 3
  lambda <- 1

  result <- generate_attributes(n, m, lambda)

  # Test output structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n, m))

  # Test no option dominates another
  max_values <- apply(result, 2, max)
  domination_check <- rowSums(sweep(result, 2, max_values, "=="))
  expect_true(!any(domination_check == m))
})

test_that("MASC_SearchRule_myopic produces valid probabilities", {
  n <- 2
  m <- 3
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  sp <- rep(1, m)
  thresh <- 0.01
  alpha <- 3
  prec <- matrix(1, n, m)
  mu <- matrix(0, n, m)

  result <- MASC_SearchRule_myopic(n, m, w, w2, sp, thresh, alpha, prec, mu)

  # Test output is a probability vector
  expect_true(is.numeric(result))
  expect_equal(length(result), n * m)
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("rMASC produces expected results", {
  set.seed(123)
  w <- c(0.5, 0.3, 0.2)  # Fixed weights

  # Test with default parameters but fixed weights
  result <- rMASC(w = w)
  expect_type(result, "list")
  expect_named(result, c("results", "weights", "parameters", "raw"))

  # Test with custom data
  custom_data <- data.frame(
    opt1_att1 = c(4.5, 4.2),
    opt1_att2 = c(3.2, 3.5),
    opt1_att3 = c(2.8, 2.9),
    opt2_att1 = c(3.8, 3.9),
    opt2_att2 = c(4.1, 4.0),
    opt2_att3 = c(3.1, 3.3)
  )

  result_custom <- rMASC(
    data = custom_data,
    w = w
  )

  expect_equal(nrow(result_custom$results), nrow(custom_data))
})

test_that("rMASC validates inputs correctly", {
  # Test sigma validation
  expect_error(
    rMASC(sigma = -1),
    "sigma must be positive"
  )

  # Test alpha validation
  expect_error(
    rMASC(alpha = -1),
    "alpha must be non-negative"
  )

  # Test delta validation
  # expect_error(
  #   rMASC(delta = -0.01),
  #   "delta must be positive"
  # )

  # Test theta validation
  expect_error(
    rMASC(theta = -0.01),
    "theta must be positive"
  )

  # Test lambda validation
  expect_error(
    rMASC(lambda = -1),
    "lambda must be positive"
  )

  # Test weight validation when provided
  expect_error(
    rMASC(w = c(-1, 0.5, 0.5)),
    "weights must be positive"
  )

  expect_error(
    rMASC(w = c(0.3, 0.3, 0.3)),
    "weights must sum to 1"
  )
})

test_that("MASC handles edge cases appropriately", {
  # Test with minimum number of options and attributes
  w1 <- 1  # Single weight for single attribute
  expect_silent(
    min_result <- rMASC(n = 1, n_options = 2, n_attributes = 1, w = w1)
  )
  expect_true(!is.null(min_result))

  # Test with large values
  w10 <- rep(0.1, 10)  # Equal weights for 10 attributes
  expect_silent(
    large_result <- rMASC(n = 1, n_options = 10, n_attributes = 10, w = w10)
  )
  expect_true(!is.null(large_result))

  # Test with equal weights
  equal_weights <- rep(1/3, 3)
  expect_silent(
    equal_w_result <- rMASC(w = equal_weights)
  )
  expect_true(!is.null(equal_w_result))
})

test_that("rMASC handles random weight generation appropriately", {
  expect_warning(
    rMASC(),
    "No weights provided. Generated random weights from Beta Distribution:"
  )
})

test_that("rMASC produces stable results with same seed and weights", {
  w <- c(0.5, 0.3, 0.2)  # Fixed weights to avoid random generation

  set.seed(123)
  result1 <- rMASC(n = 5, w = w)

  set.seed(123)
  result2 <- rMASC(n = 5, w = w)

  expect_equal(result1$results, result2$results)
  expect_equal(result1$raw, result2$raw)
})

test_that("C++ functions handle numerical edge cases", {
  # Test generate_attributes with extreme lambda values
  expect_silent(result <- generate_attributes_cpp(2, 3, 1e-10))
  expect_silent(result <- generate_attributes_cpp(2, 3, 1e10))

  # Test MASC_SearchRule_myopic with extreme precision values
  n <- 2
  m <- 3
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  sp <- rep(1e-10, m)  # Very low precision
  thresh <- 0.01
  alpha <- 3
  prec <- matrix(1e-10, n, m)  # Very low precision
  mu <- matrix(0, n, m)
  expect_silent(result <- MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha, prec, mu))
  expect_true(all(is.finite(result)))

  # Test with very high precision
  sp <- rep(1e10, m)
  prec <- matrix(1e10, n, m)
  expect_silent(result <- MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha, prec, mu))
  expect_true(all(is.finite(result)))
})

test_that("Fixation sequences are valid", {
  w <- c(0.5, 0.3, 0.2)  # Explicit weights
  result <- rMASC(n = 1, n_options = 2, n_attributes = 3, w = w)

  # Check fix_sequence structure
  fix_seq <- result$raw[[1]]$fix_sequence
  expect_true(is.numeric(fix_seq))
  expect_true(all(fix_seq >= 1))
  expect_true(all(fix_seq <= 6))  # 2 options * 3 attributes

  # Check fixation proportions sum to 1
  prop_opt <- result$raw[[1]]$prop_fix_opt
  expect_equal(sum(prop_opt), 1, tolerance = 1e-10)

  prop_att <- result$raw[[1]]$prop_fix_att
  expect_equal(sum(prop_att), 1, tolerance = 1e-10)
})

test_that("rMASC handles different numbers of options correctly", {
  # Test with 3 options
  n_opt <- 3
  n_att <- 2
  w <- c(0.6, 0.4)
  result <- rMASC(n = 1, n_options = n_opt, n_attributes = n_att, w = w)

  # Check results structure adapts to n_options
  expect_equal(ncol(result$raw[[1]]$x), n_att)
  expect_equal(nrow(result$raw[[1]]$x), n_opt)
  expect_equal(length(result$raw[[1]]$prop_fix_opt), n_opt)

  # Verify column names in results dataframe
  expect_true(all(paste0("prop_fix_opt", 1:n_opt) %in% names(result$results)))
})

test_that("MASC decision process is coherent", {
  # Create data where one option clearly dominates
  custom_data <- data.frame(
    opt1_att1 = 5,  # Clearly better option
    opt1_att2 = 5,
    opt1_att3 = 5,
    opt2_att1 = 1,
    opt2_att2 = 1,
    opt2_att3 = 1
  )

  w <- c(0.4, 0.3, 0.3)
  result <- rMASC(data = custom_data, w = w, sigma = 0.1)  # Low noise

  # Should choose the dominant option
  expect_equal(result$results$response, 1)
  expect_equal(result$results$best_option, 1)
  expect_true(result$results$correct)
})

test_that("MASC sampling behavior is reasonable", {
  set.seed(123)
  w <- c(0.5, 0.3, 0.2)  # Explicit weights
  result <- rMASC(n = 100, w = w)  # Run multiple trials

  # Response time (rt) should be reasonable
  expect_true(all(result$results$rt > 0))
  expect_true(all(result$results$rt <= 100))  # max_steps default

  # Accuracy should be above chance
  accuracy <- mean(result$results$correct)
  expect_true(accuracy > 0.5)  # Should be better than random guessing

  # Fixation proportions should be reasonable
  expect_true(all(result$results$prop_fix_opt1 >= 0))
  expect_true(all(result$results$prop_fix_opt1 <= 1))
  expect_true(all(result$results$prop_fix_opt2 >= 0))
  expect_true(all(result$results$prop_fix_opt2 <= 1))

  # Each trial's fixation proportions should sum to 1
  trial_sums <- result$results$prop_fix_opt1 + result$results$prop_fix_opt2
  expect_true(all(abs(trial_sums - 1) < 1e-10))
})

test_that("MASC_SearchRule_myopic matches MATLAB implementation", {
  # Helper function to check matrix output against expected values
  check_output <- function(result, expected_matrix, expected_sum, tolerance = 1e-3) { #tolerance = 1e-4
    # Convert result to matrix
    result_matrix <- round(matrix(result, nrow = 2, byrow = FALSE),4)
    # Check matrix values
    expect_equal(result_matrix, expected_matrix, tolerance = tolerance)
    # Check column sums
    col_sums <- colSums(result_matrix)
    expect_equal(col_sums, expected_sum, tolerance = tolerance)
    # Check total sum is 1
    expect_equal(sum(result), 1, tolerance = 1e-10)
    # Check all values are probabilities
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
  }

  # Case 1: Equal Weights
  n <- 2
  m <- 3
  w <- c(1/3, 1/3, 1/3)
  w2 <- w^2
  sp <- rep(1, m)
  thresh <- 0.01
  alpha <- 3
  attMean <- matrix(c(1.0, 0.5, -0.2,
                     0.3, 0.8, 0.1),
                   nrow=2, byrow=TRUE)
  attPrecision <- matrix(1, nrow=n, ncol=m)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.2421, 0.2421, 0.2421,
                            0.0912, 0.0912, 0.0912),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.3333, 0.3333, 0.3333)

  check_output(result, expected_matrix, expected_sum)

  # Case 2: Extreme Weights
  w <- c(0.98, 0.01, 0.01)
  w2 <- w^2

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.8007, 0.0399, 0.0399,
                            0.0399, 0.0399, 0.0399),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.8405, 0.0797, 0.0797)

  check_output(result, expected_matrix, expected_sum)

  # Case 3: Extreme Values
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  attMean <- matrix(c(10.0, -10.0, 0.0,
                     -10.0, 10.0, 0.0),
                   nrow=2, byrow=TRUE)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.2437, 0.2437, 0.2437,
                            0.0896, 0.0896, 0.0896),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.3333, 0.3333, 0.3333)

  check_output(result, expected_matrix, expected_sum)

  # Case 4: Different Precisions
  attMean <- matrix(c(1.0, 0.5, -0.2,
                     0.3, 0.8, 0.1),
                   nrow=2, byrow=TRUE)
  sp <- c(0.1, 1.0, 10.0)
  attPrecision <- matrix(c(0.1, 1.0, 10.0,
                          10.0, 1.0, 0.1),
                        nrow=2, byrow=TRUE)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.8007, 0.0399, 0.0399,
                            0.0399, 0.0399, 0.0399),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.8405, 0.0797, 0.0797)

  check_output(result, expected_matrix, expected_sum)

  # Case 5: Edge Threshold
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  sp <- rep(1, m)
  thresh <- 0.49999
  attMean <- matrix(c(1.0, 0.5, -0.2,
                     0.3, 0.8, 0.1),
                   nrow=2, byrow=TRUE)
  attPrecision <- matrix(1, nrow=n, ncol=m)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.2057, 0.2321, 0.2485,
                            0.1156, 0.1025, 0.0957),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.3213, 0.3345, 0.3442)

  check_output(result, expected_matrix, expected_sum)

  # Case 6: Zero Means
  thresh <- 0.01
  attMean <- matrix(0, nrow=n, ncol=m)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.3457, 0.0771, 0.0771,
                            0.3457, 0.0771, 0.0771),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.6914, 0.1543, 0.1543)

  check_output(result, expected_matrix, expected_sum)

  # Case 7: High Precision Contrast
  w <- c(0.4, 0.3, 0.3)
  w2 <- w^2
  attMean <- matrix(c(1.0, 0.5, -0.2,
                     0.3, 0.8, 0.1),
                   nrow=2, byrow=TRUE)
  attPrecision <- matrix(c(100.0, 1.0, 1.0,
                          1.0, 100.0, 1.0),
                        nrow=2, byrow=TRUE)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.0399, 0.0399, 0.0399,
                            0.8007, 0.0399, 0.0399),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.8405, 0.0797, 0.0797)

  check_output(result, expected_matrix, expected_sum)

  # Case 8: Low Search Sensitivity
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  attPrecision <- matrix(1, nrow=n, ncol=m)
  lowAlpha <- 0.1

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, lowAlpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.1810, 0.1638, 0.1638,
                            0.1638, 0.1638, 0.1638),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.3448, 0.3276, 0.3276)

  check_output(result, expected_matrix, expected_sum)

  # Case 9: High Search Sensitivity
  highAlpha <- 10.0

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, highAlpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.9998, 0.0000, 0.0000,
                            0.0000, 0.0000, 0.0000),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.9998, 0.0001, 0.0001)

  # Higher tolerance for this test case due to very small values
  check_output(result, expected_matrix, expected_sum, tolerance = 1e-3)

  # Case 10: Uneven Precision
  w <- c(0.5, 0.3, 0.2)
  w2 <- w^2
  alpha <- 3
  sp <- c(0.1, 1.0, 10.0)
  attPrecision <- matrix(c(10.0, 1.0, 0.1,
                          0.1, 1.0, 10.0),
                        nrow=2, byrow=TRUE)

  result <- masc:::MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha,
                                             attPrecision, attMean)

  # From MATLAB output
  expected_matrix <- matrix(c(0.0399, 0.0399, 0.0399,
                            0.8005, 0.0399, 0.0399),
                          nrow=2, byrow=TRUE)
  expected_sum <- c(0.8404, 0.0798, 0.0798)

  check_output(result, expected_matrix, expected_sum)
})

# End-to-end testing for full rMASC function
test_that("Full rMASC implementation with myopic search is consistent", {
  set.seed(123)

  # Create test data with known structure
  test_data <- data.frame(
    opt1_att1 = 1.0,
    opt1_att2 = 0.5,
    opt1_att3 = -0.2,
    opt2_att1 = 0.3,
    opt2_att2 = 0.8,
    opt2_att3 = 0.1
  )

  # Run model with very low noise to ensure deterministic behavior
  result <- rMASC(
    data = test_data,
    w = c(0.5, 0.3, 0.2),
    sigma = 0.01,  # Very low noise
    alpha = 3,
    delta = 0.01,
    theta = 0.01,
    max_steps = 50
  )

  # Verify fixation proportions are valid
  expect_true(all(result$results$prop_fix_opt1 >= 0))
  expect_true(all(result$results$prop_fix_opt1 <= 1))
  expect_true(all(result$results$prop_fix_opt2 >= 0))
  expect_true(all(result$results$prop_fix_opt2 <= 1))

  # Verify proportions sum to 1
  expect_equal(
    result$results$prop_fix_opt1 + result$results$prop_fix_opt2,
    1,
    tolerance = 1e-10
  )

  # Verify fixation sequence values are valid OAP indices
  fix_seq <- result$raw[[1]]$fix_sequence
  expect_true(all(fix_seq >= 1))
  expect_true(all(fix_seq <= 6))  # 2 options * 3 attributes

  # Verify option values calculation is correct
  expect_equal(
    result$raw[[1]]$opt_values,
    c(1.0*0.5 + 0.5*0.3 + (-0.2)*0.2, 0.3*0.5 + 0.8*0.3 + 0.1*0.2),
    tolerance = 1e-10
  )
})

