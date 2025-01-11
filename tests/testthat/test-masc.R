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
  expect_error(
    rMASC(delta = -0.01),
    "delta must be positive"
  )

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


