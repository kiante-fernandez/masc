#' Generate Non-dominated Attribute Values
#'
#' @description
#' Internal function that generates attribute values for options ensuring no option
#' dominates another (i.e., is better on all attributes).
#'
#' @param n Integer. Number of options.
#' @param m Integer. Number of attributes per option.
#' @param lambda Numeric. Precision parameter for attribute value generation.
#'
#' @return Matrix of attribute values (n x m).
#' @keywords internal
generate_attributes <- function(n, m, lambda) {
  generate_attributes_cpp(n, m, lambda)
}

#' Myopic Search Rule for MASC Model
#'
#' @description
#' Internal function implementing the myopic search rule for the MASC model.
#' Calculates transition probabilities for the next fixation by evaluating how
#' likely each possible sample would lead to a decision.
#'
#' @param n Integer. Number of options.
#' @param m Integer. Number of attributes.
#' @param w Numeric vector. Attribute weights.
#' @param w2 Numeric vector. Squared attribute weights (pre-computed).
#' @param sp Numeric vector. Sampling precision for each attribute.
#' @param thresh Numeric. Current decision threshold.
#' @param alpha Numeric. Search sensitivity parameter.
#' @param prec Matrix. Current attribute precision (n x m).
#' @param mu Matrix. Current attribute means (n x m).
#'
#' @return Numeric vector. Probabilities for each possible fixation location.
#' @keywords internal
MASC_SearchRule_myopic <- function(n, m, w, w2, sp, thresh, alpha, prec, mu) {
  # Call C++ implementation
  result <- MASC_SearchRule_myopic_cpp(n, m, w, w2, sp, thresh, alpha, prec, mu)
  return(result)
}

#' Multi-Attribute Search and Choice (MASC) Model
#'
#' @description
#' Implements the MASC model of multi-attribute decision making. This model simulates
#' how people make decisions when comparing options with multiple attributes by
#' sequentially sampling information about different attributes until reaching a decision.
#'
#' @param data Optional data frame containing trial-wise attribute values. Each row
#'   represents one trial, and columns should be named following the pattern
#'   "opt{i}_att{j}" where i is the option number and j is the attribute number.
#'   For example, with 2 options and 3 attributes, columns should be:
#'   opt1_att1, opt1_att2, opt1_att3, opt2_att1, opt2_att2, opt2_att3.
#'   If NULL, generates random values for n trials.
#' @param n Integer. Number of trials to generate when data is NULL (default: 1).
#'   Ignored when data is provided.
#' @param n_options Integer. Number of choice options (default: 2).
#' @param n_attributes Integer. Number of attributes per option (default: 3).
#' @param w Numeric vector. Attribute weights summing to 1. If NULL, weights are
#'   randomly generated from beta(3/4, 3/4) distribution (default: NULL).
#' @param sigma Numeric. Standard deviation of sampling noise (default: 1).
#' @param alpha Numeric. Controls how strongly fixations follow the myopic search
#'   rule. Higher values (>10) make search more deterministic, lower values (≈0)
#'   make it more random (default: 3).
#' @param delta Numeric. Amount by which decision threshold increases per fixation
#'   (default: 0.01).
#' @param theta Numeric. Initial decision threshold (default: 0.01).
#' @param lambda Numeric. Precision of prior beliefs about attributes (default: 1).
#' @param max_steps Integer. Maximum number of fixations allowed (default: 100).
#'
#' @return A list containing:
#' \itemize{
#'   \item results: Data frame with trial-by-trial results including:
#'     \itemize{
#'       \item trial: Trial number
#'       \item response: Option chosen by model (1 to n_options)
#'       \item best_option: Option with highest weighted value
#'       \item correct: Whether response matches best_option
#'       \item rt: Number of fixations taken
#'       \item value_difference: Difference between best and worst option values
#'       \item prop_fix_opt1, prop_fix_opt2: Proportion of fixations to each option
#'     }
#'   \item weights: Vector of attribute weights used
#'   \item parameters: List of model parameters used (sigma, alpha, delta, theta)
#'   \item raw: List containing detailed raw data for each trial
#' }
#'
#' @examples
#' # Example 1: Generate 5 random trials
#' results <- rMASC(n = 5)
#'
#' # Example 2: Custom attribute values for multiple trials
#' trial_data <- data.frame(
#'     # Option 1's attributes across 3 trials
#'     opt1_att1 = c(4.5, 4.2, 4.8),  # Attribute 1 values
#'     opt1_att2 = c(3.2, 3.5, 3.1),  # Attribute 2 values
#'     opt1_att3 = c(2.8, 2.9, 2.7),  # Attribute 3 values
#'     # Option 2's attributes across 3 trials
#'     opt2_att1 = c(3.8, 3.9, 3.7),  # Attribute 1 values
#'     opt2_att2 = c(4.1, 4.0, 4.2),  # Attribute 2 values
#'     opt2_att3 = c(3.1, 3.3, 3.0)   # Attribute 3 values
#' )
#'
#' # Run model with custom weights
#' results <- rMASC(
#'     data = trial_data,
#'     w = c(0.5, 0.3, 0.2)  # weights for attributes
#' )
#'
#' @references
#' Gluth, S., Deakin, J., & Rieskamp, J. (2024). A Theory of Multi-Attribute Search
#' and Choice.
#'
#' @export
rMASC <- function(data = NULL,
                  n = 1,
                  n_options = 2,
                  n_attributes = 3,
                  w = NULL,
                  sigma = 1,
                  alpha = 3,
                  delta = 0.01,
                  theta = 0.01,
                  lambda = 1,
                  max_steps = 100) {

  # Validate numeric parameters
  if (sigma <= 0) stop("sigma must be positive")
  if (alpha < 0) stop("alpha must be non-negative")
  if (delta <= 0) stop("delta must be positive")
  if (theta <= 0) stop("theta must be positive")
  if (lambda <= 0) stop("lambda must be positive")
  if(!is.numeric(n) || n < 1 || n != round(n))
    stop("n must be a positive integer")

  # Determine number of trials
  n_trials <- if(!is.null(data)) {
    if(!is.data.frame(data)) stop("data must be a data frame")
    nrow(data)
  } else {
    n
  }

  # Generate or validate weights
  if(is.null(w)) {
    w <- rbeta(n_attributes, 3/4, 3/4)
    w <- w/sum(w)
    warning("No weights provided. Generated random weights from Beta Distribution: ",
            paste(round(w, 3), collapse=", "))
  }

  # Validate weights
  if(length(w) != n_attributes)
    stop("Length of weights must match number of attributes")
  if(any(w <= 0))
    stop("weights must be positive")
  if(abs(sum(w) - 1) > .Machine$double.eps)
    stop("weights must sum to 1")

  # If data provided, validate column names
  if(!is.null(data)) {
    expected_cols <- c()
    for(i in 1:n_options) {
      for(j in 1:n_attributes) {
        expected_cols <- c(expected_cols, sprintf("opt%d_att%d", i, j))
      }
    }
    missing_cols <- setdiff(expected_cols, names(data))
    if(length(missing_cols) > 0) {
      stop("Missing columns: ", paste(missing_cols, collapse=", "),
           "\nColumn names should be opt1_att1, opt1_att2, etc.")
    }
  }

  # Pre-allocate results
  all_trials <- vector("list", n_trials)

  # Process each trial
  for(trial in 1:n_trials) {
    # Get or generate stimulus values
    trial_x <- if(is.null(data)) {
      generate_attributes_cpp(n_options, n_attributes, lambda)
    } else {
      matrix(as.numeric(data[trial, ]), nrow=n_options, byrow=TRUE)
    }

    # Run sampling process using C++
    trial_results <- rMASC_sampling_cpp(
      trial_x = trial_x,
      w = w,
      sigma = sigma,
      alpha = alpha,
      delta = delta,
      theta = theta,
      lambda = lambda,
      max_steps = max_steps,
      n_options = n_options,
      n_attributes = n_attributes
    )

    # Calculate option values
    opt_values <- drop(trial_x %*% w)

    # Store trial results
    all_trials[[trial]] <- list(
      trial = trial,
      response = which.max(trial_results$response),
      best_option = trial_results$best_option,
      correct = which.max(trial_results$response) == trial_results$best_option,
      rt = trial_results$rt,
      x = trial_x,
      opt_values = opt_values,
      weights = w,
      sigma = sigma,
      alpha = alpha,
      delta = delta,
      theta = theta,
      fix_sequence = trial_results$fix_sequence,
      prop_fix_opt = trial_results$prop_fix_opt,
      prop_fix_att = trial_results$prop_fix_att
    )
  }

  # Convert results to tidy data frame
  results_df <- do.call(rbind, lapply(all_trials, function(trial) {
    # Create base results
    base_results <- data.frame(
      trial = trial$trial,
      response = trial$response,
      best_option = trial$best_option,
      correct = trial$correct,
      rt = trial$rt
    )

    # Add proportion fixation columns for each option
    for(i in 1:n_options) {
      base_results[[paste0("prop_fix_opt", i)]] <- trial$prop_fix_opt[i]
    }

    base_results
  }))

  # Add original trial data if provided
  if(!is.null(data)) {
    results_df <- cbind(results_df, data)
  }

  # Return results
  return(list(
    results = results_df,  # Main results in tidy format
    weights = w,           # Weights used
    parameters = list(     # Model parameters used
      sigma = sigma,
      alpha = alpha,
      delta = delta,
      theta = theta
    ),
    raw = all_trials       # Raw trial data if needed
  ))
}

