
<!-- README.md is generated from README.Rmd. Please edit that file -->

# masc: Multi-Attribute Search and Choice Model in R

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

The `masc` package implements the Multi-Attribute Search and Choice
(MASC) model, a hierarchical Bayesian framework for understanding how
people make decisions between options with multiple attributes. Based on
the work by Gluth, Deakin, & Rieskamp (2024), this package simulates:

- Information search patterns in multi-attribute decisions
- Belief updating about attribute values
- Choice dynamics and decision termination
- The interplay between attention and valuation

## Installation

``` r
# Install development version from GitHub
devtools::install_github("kiantefernandez/masc")
```

## Basic Usage

### Simple Random Trials

``` r
library(masc)

# Generate 5 random trials
results <- rMASC(
  n = 5,
  w = c(0.5, 0.3, 0.2)  # weights for attributes
)

# View results
print(results$results)
```

### Custom Attribute Values

``` r
# Create trial data
trial_data <- data.frame(
  # Option 1's attributes across 3 trials
  opt1_att1 = c(4.5, 4.2, 4.8),  # Attribute 1 values
  opt1_att2 = c(3.2, 3.5, 3.1),  # Attribute 2 values
  opt1_att3 = c(2.8, 2.9, 2.7),  # Attribute 3 values
  # Option 2's attributes across 3 trials
  opt2_att1 = c(3.8, 3.9, 3.7),  # Attribute 1 values
  opt2_att2 = c(4.1, 4.0, 4.2),  # Attribute 2 values
  opt2_att3 = c(3.1, 3.3, 3.0)   # Attribute 3 values
)

# Run model with custom data
results <- rMASC(
  data = trial_data,
  w = c(0.5, 0.3, 0.2)  # weights for attributes
)
```

## Function Parameters

The `rMASC()` function accepts the following parameters:

- `data`: Optional data frame with trial-wise attribute values
- `n`: Number of trials to generate when `data` is NULL (default: 1)
- `n_options`: Number of choice options (default: 2)
- `n_attributes`: Number of attributes per option (default: 3)
- `w`: Attribute weights summing to 1 (if NULL, randomly generated)
- `sigma`: Standard deviation of sampling noise (default: 1)
- `alpha`: Search rule sensitivity (default: 3)
- `delta`: Decision threshold increase per fixation (default: 0.01)
- `theta`: Initial decision threshold (default: 0.01)
- `lambda`: Precision of prior beliefs about attributes (default: 1)
- `max_steps`: Maximum number of fixations allowed (default: 100)

## Output Structure

The function returns a list containing:

- `results`: Data frame with trial-by-trial results including:
  - trial number
  - chosen option
  - best option
  - choice consistency
  - number of fixations (rt)
  - fixation proportions
- `weights`: Vector of attribute weights used
- `parameters`: List of model parameters used
- `raw`: List containing detailed raw data for each trial

The `raw` component contains a list where each element corresponds to a
trial and includes:

- `trial`: Trial number
- `response`: The option chosen by the model (1 to n_options)
- `best_option`: The option with the highest weighted value
- `correct`: Boolean indicating if response matches best_option
- `rt`: Number of fixations taken to reach a decision
- `x`: Matrix of true attribute values for all options
- `opt_values`: Vector of computed option values (weighted sums)
- `weights`: Vector of attribute weights used in this trial
- `sigma`: Sampling noise parameter used
- `alpha`: Search sensitivity parameter used
- `delta`: Threshold increment parameter used
- `theta`: Initial threshold parameter used
- `fix_sequence`: Vector showing the sequence of fixations made
- `prop_fix_opt`: Vector of proportions of fixations to each option
- `prop_fix_att`: Vector of proportions of fixations to each attribute

## Reference

Gluth, S., Deakin, J., & Rieskamp, J. (2024). A Theory of
Multi-Attribute Search and Choice.
<https://doi.org/10.31234/osf.io/3qzak>

## License

This package is released under the MIT License.
