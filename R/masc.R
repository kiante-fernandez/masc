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

# ---------------------------------------------------------------------------
# Internal helpers for correlation/covariance structures (not exported)
# ---------------------------------------------------------------------------

# Coerce a user-supplied structure to a valid m x m matrix:
#   NULL            -> identity (when default_identity) else NULL
#   single number 0 -> identity (independent attributes)
#   single number r -> equicorrelation matrix (1 on diagonal, r off-diagonal)
#   matrix          -> used as-is (size-validated)
# The result is symmetrized and made positive-definite via .masc_ensure_pd().
.masc_resolve_sigma <- function(Sigma, m, default_identity = TRUE) {
  if (is.null(Sigma)) {
    if (default_identity) return(diag(m)) else return(NULL)
  }
  if (!is.matrix(Sigma)) {
    if (length(Sigma) != 1 || !is.numeric(Sigma))
      stop("Sigma must be a single number or an ", m, " x ", m, " matrix")
    rho <- Sigma
    if (rho == 0) return(diag(m))
    Sigma <- matrix(rho, m, m)
    diag(Sigma) <- 1
  }
  if (nrow(Sigma) != m || ncol(Sigma) != m)
    stop("Sigma must be an ", m, " x ", m, " matrix")
  .masc_ensure_pd(Sigma)
}

# Symmetrize and clamp eigenvalues to keep the matrix positive-definite.
.masc_ensure_pd <- function(Sigma, eps = 1e-10) {
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  if (any(eig$values < eps)) {
    eig$values[eig$values < eps] <- eps
    Sigma <- eig$vectors %*% diag(eig$values, length(eig$values)) %*% t(eig$vectors)
    Sigma <- (Sigma + t(Sigma)) / 2
  }
  Sigma
}

.masc_is_identity <- function(M, tol = 1e-10) {
  is.matrix(M) && nrow(M) == ncol(M) && max(abs(M - diag(nrow(M)))) < tol
}

.masc_is_diagonal <- function(M, tol = 1e-10) {
  is.matrix(M) && max(abs(M - diag(diag(M)))) < tol
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
#'   "opt<i>_att<j>" where i is the option number and j is the attribute number.
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
#'   rule. Higher values (>10) make search more deterministic, lower values
#'   (near 0) make it more random (default: 3).
#' @param delta Numeric. Amount by which decision threshold increases per fixation
#'   (default: 0.01).
#' @param theta Numeric. Initial decision threshold (default: 0.01).
#' @param lambda Numeric. Precision of prior beliefs about attributes (default: 1).
#' @param max_steps Integer. Maximum number of fixations allowed (default: 100).
#' @param Sigma_true Correlation/covariance structure of the generated stimuli
#'   (an `n_attributes` x `n_attributes` matrix, or a single number giving a
#'   uniform off-diagonal correlation). When `NULL` (default), attributes are
#'   independent (identity), reproducing the original MASC behaviour. Ignored
#'   when `data` is supplied (the data are the stimuli).
#' @param Sigma_belief The decision maker's assumed correlation structure between
#'   attributes (matrix or single number). This is what enables the multivariate
#'   ("MASC-C") belief update: observing one attribute spreads information to
#'   correlated attributes via a Kalman update. `NULL` (default) matches
#'   `Sigma_true`; `0` forces independent (univariate) beliefs. When the
#'   resulting matrix is diagonal the model reduces exactly to the original
#'   univariate MASC update.
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
#'       \item prop_fix_opt1, prop_fix_opt2: Proportion of fixations to each option
#'     }
#'   \item weights: Vector of attribute weights used
#'   \item parameters: List of model parameters used (sigma, alpha, delta, theta)
#'   \item raw: List containing detailed raw data for each trial. Each element corresponds to a trial and includes:
#'     \itemize{
#'       \item trial: Trial number
#'       \item response: The option chosen by the model (1 to n_options)
#'       \item best_option: The option with the highest weighted value
#'       \item correct: Boolean indicating if response matches best_option
#'       \item rt: Number of fixations taken to reach a decision
#'       \item x: Matrix of true attribute values for all options
#'       \item opt_values: Vector of computed option values (weighted sums)
#'       \item weights: Vector of attribute weights used in this trial
#'       \item sigma: Sampling noise parameter used
#'       \item alpha: Search sensitivity parameter used
#'       \item delta: Threshold increment parameter used
#'       \item theta: Initial threshold parameter used
#'       \item fix_sequence: Vector showing the sequence of fixations made
#'       \item prop_fix_opt: Vector of proportions of fixations to each option
#'       \item prop_fix_att: Vector of proportions of fixations to each attribute
#'     }
#' }
#'
#' @examples
#' # Example 1: Generate 5 random trials
#' results <- rMASC(n = 5, w = c(0.5, 0.3, 0.2))
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
#' # Example 3: Correlated attributes (MASC-C). The decision maker exploits a
#' # positive correlation structure, so observing one attribute informs beliefs
#' # about the others ("belief spread").
#' results <- rMASC(
#'     n = 20,
#'     w = c(0.5, 0.3, 0.2),
#'     Sigma_true = 0.6,    # stimuli are positively correlated
#'     Sigma_belief = 0.6   # matched beliefs (use 0 for the original MASC model)
#' )
#'
#' @references
#' Gluth, S., Deakin, J., & Rieskamp, J. (2026). A theory of multiattribute search
#' and choice. *Psychological Review*. <https://doi.org/10.1037/rev0000614>
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
                  max_steps = 100,
                  Sigma_true = NULL,
                  Sigma_belief = NULL) {
  # Validate numeric parameters
  if (sigma <= 0) stop("sigma must be positive")
  if (alpha < 0) stop("alpha must be non-negative")
#  if (delta <= 0) stop("delta must be positive")
  if (theta <= 0) stop("theta must be positive")
  if (lambda <= 0) stop("lambda must be positive")
  if(!is.numeric(n) || n < 1 || n != round(n))
    stop("n must be a positive integer")

  # Convert scalar sigma to vector if needed
  if(length(sigma) == 1) {
    sigma <- rep(sigma, n_attributes)
  } else if(length(sigma) != n_attributes) {
    stop("Length of sigma vector must match number of attributes")
  }

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

  # Resolve correlation structures (see .masc_resolve_sigma / .masc_ensure_pd).
  # Sigma_true governs the generated stimuli; Sigma_belief governs the (possibly
  # multivariate) belief update. NULL/identity/diagonal reduce to original MASC.
  Sigma_true <- .masc_resolve_sigma(Sigma_true, n_attributes, default_identity = TRUE)
  if(is.null(Sigma_belief)) {
    Sigma_belief <- Sigma_true            # matched beliefs
  } else {
    Sigma_belief <- .masc_resolve_sigma(Sigma_belief, n_attributes, default_identity = TRUE)
  }
  # Only pass non-trivial structures down to C++ so the default path is byte-identical:
  # an identity Sigma_true uses the original stimulus draw; a diagonal Sigma_belief
  # uses the original univariate update.
  Sigma_true_arg   <- if(.masc_is_identity(Sigma_true)) NULL else Sigma_true
  Sigma_belief_arg <- if(.masc_is_diagonal(Sigma_belief)) NULL else Sigma_belief

  # If data provided, validate column names
  if(!is.null(data)) {
    expected_cols <- outer(1:n_options, 1:n_attributes,
                           FUN = function(i, j) sprintf("opt%d_att%d", i, j))
    expected_cols <- as.vector(expected_cols)
    missing_cols <- setdiff(expected_cols, names(data))
    if(length(missing_cols) > 0) {
      stop("Missing columns: ", paste(missing_cols, collapse=", "),
           "\nColumn names should be opt1_att1, opt1_att2, etc.")
    }
  }

  # Pre-compute column names for fixation proportions
  fix_opt_cols <- paste0("prop_fix_opt", seq_len(n_options))

  # Pre-allocate results data frame with all necessary columns
  results_df <- data.frame(
    trial = seq_len(n_trials),
    response = integer(n_trials),
    best_option = integer(n_trials),
    correct = logical(n_trials),
    rt = integer(n_trials)
  )

  # Pre-allocate fixation proportion matrix
  prop_fix_matrix <- matrix(0, nrow = n_trials, ncol = n_options,
                            dimnames = list(NULL, fix_opt_cols))

  # Pre-allocate raw results list
  all_trials <- vector("list", n_trials)

  # Pre-allocate trial data if not provided
  if(is.null(data)) {
    trial_data <- replicate(n_trials,
                            generate_attributes_cpp(n_options, n_attributes, lambda,
                                                    Sigma_true_arg),
                            simplify = FALSE)
  }

  # Process each trial
  for(trial in seq_len(n_trials)) {
    # Get or extract stimulus values
    trial_x <- if(is.null(data)) {
      trial_data[[trial]]
    } else {
      # Check dimensions before creating matrix
      if (length(data[trial, ]) != n_options * n_attributes) {
        stop(sprintf("Error in trial %d: Number of data values (%d) does not match expected values for %d options and %d attributes (%d)",
                     trial, length(data[trial, ]), n_options, n_attributes, n_options * n_attributes))
      }
      matrix(as.numeric(data[trial, ]), nrow=n_options, byrow=TRUE)
    }

    # After creating matrix, check dimensions explicitly
    if (nrow(trial_x) != n_options || ncol(trial_x) != n_attributes) {
      stop(sprintf("Error in trial %d: Matrix dimensions (%d, %d) do not match expected dimensions (%d, %d)",
                   trial, nrow(trial_x), ncol(trial_x), n_options, n_attributes))
    }
    #print(dim(trial_x))
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
      n_attributes = n_attributes,
      Sigma_belief = Sigma_belief_arg
    )

    # Calculate option values
    opt_values <- drop(trial_x %*% w)
    # Efficient batch assignment of trial results
    results_df$response[trial] <- which.max(trial_results$response)
    results_df$best_option[trial] <- trial_results$best_option
    results_df$rt[trial] <- trial_results$rt
    results_df$correct[trial] <-
      which.max(trial_results$response) == trial_results$best_option

    # Store fixation proportions in matrix
    prop_fix_matrix[trial, ] <- trial_results$prop_fix_opt

    # Store full results in raw list
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

  # Combine results_df with prop_fix_matrix
  results_df <- cbind(results_df, prop_fix_matrix)

  # Add original trial data if provided
  if(!is.null(data)) {
    results_df <- cbind(results_df, data)
  }

  # Return results
  return(list(
    results = results_df,
    weights = w,
    parameters = list(
      sigma = sigma,
      alpha = alpha,
      delta = delta,
      theta = theta,
      Sigma_true = Sigma_true,
      Sigma_belief = Sigma_belief
    ),
    raw = all_trials
  ))
}

