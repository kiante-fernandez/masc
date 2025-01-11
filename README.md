
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

The primary function `rMASC()` allows you to simulate the MASC model:

### Basic Example: Random Trials

``` r
library(masc)

# Generate 5 random trials
results <- rMASC(n = 5)
print(results$results)
```

### Custom Attribute Values

``` r
# Create a data frame with custom attribute values
trial_data <- data.frame(
  opt1_att1 = c(4.5, 4.2, 4.8),  # Option 1, Attribute 1 values
  opt1_att2 = c(3.2, 3.5, 3.1),  # Option 1, Attribute 2 values
  opt1_att3 = c(2.8, 2.9, 2.7),  # Option 1, Attribute 3 values
  opt2_att1 = c(3.8, 3.9, 3.7),  # Option 2, Attribute 1 values
  opt2_att2 = c(4.1, 4.0, 4.2),  # Option 2, Attribute 2 values
  opt2_att3 = c(3.1, 3.3, 3.0)   # Option 2, Attribute 3 values
)

# Run model with custom attribute values and weights
results <- rMASC(
  data = trial_data,
  w = c(0.5, 0.3, 0.2)  # custom attribute weights
)
```

### Simple Decision Simulation

``` r
library(masc)

# Simulate 100 decisions between two options with three attributes each
results <- rMASC(
  n = 100,
  w = c(0.5, 0.3, 0.2)  # Attribute weights
)

# Examine basic results
head(results$results)
```

### Custom Decision Scenarios

``` r
# Create a choice scenario between two smartphones
phone_data <- data.frame(
  # Option 1: High-end phone
  opt1_att1 = c(4.5),  # Screen size (inches)
  opt1_att2 = c(5000), # Battery (mAh)
  opt1_att3 = c(256),  # Storage (GB)
  
  # Option 2: Mid-range phone
  opt2_att1 = c(4.2),  
  opt2_att2 = c(4000),
  opt2_att3 = c(128)
)

# Run model with standardized values and specific parameters
phone_choice <- rMASC(
  data = scale(phone_data),  # Standardize attributes
  w = c(0.4, 0.4, 0.2),     # Equal weight to screen and battery
  alpha = 5,                 # More deterministic search
  sigma = 0.5               # Lower sampling noise
)
```

### Analyzing Search Patterns

``` r
# Simulate decisions with different attribute importance
high_dispersion <- rMASC(
  n = 50,
  w = c(0.8, 0.1, 0.1),  # One dominant attribute
  alpha = 3
)

low_dispersion <- rMASC(
  n = 50,
  w = c(0.4, 0.3, 0.3),  # More equal weights
  alpha = 3
)

# Compare fixation patterns
mean(high_dispersion$results$rt)     # Should be lower
mean(low_dispersion$results$rt)      # Should be higher
```

## Key Features

The MASC model captures several key empirical phenomena:

1.  **Choice Consistency**: More consistent choices when value
    differences are large

``` r
# Value difference affects choice consistency
easy_choice <- rMASC(
  data = data.frame(
    opt1_att1 = 2, opt1_att2 = 2, opt1_att3 = 2,
    opt2_att1 = 0, opt2_att2 = 0, opt2_att3 = 0
  ),
  w = c(0.4, 0.3, 0.3)
)
```

2.  **Search Efficiency**: Attention allocation follows attribute
    importance

``` r
# More important attributes are sampled first
results <- rMASC(
  n = 1,
  w = c(0.7, 0.2, 0.1),
  alpha = 5
)
```

3.  **Attraction Search Effect**: Tendency to keep sampling from
    promising options

``` r
# Higher alpha increases attraction search effect
results <- rMASC(
  n = 1,
  alpha = 8,  # Strong attraction to promising options
  w = c(0.4, 0.3, 0.3)
)
```

## Function Parameters

The `rMASC()` function accepts the following parameters:

- `data`: Optional data frame with trial-wise attribute values
- `n`: Number of trials to generate when `data` is NULL (default: 1)
- `n_options`: Number of choice options (default: 2)
- `n_attributes`: Number of attributes per option (default: 3)
- `w`: Attribute weights (if NULL, randomly generated)
- `sigma`: Standard deviation of sampling noise (default: 1)
- `alpha`: Search sensitivity parameter (default: 3)
- `delta`: Decision threshold increase per fixation (default: 0.01)
- `theta`: Initial decision threshold (default: 0.01)
- `lambda`: Precision of prior beliefs about attributes (default: 1)
- `max_steps`: Maximum number of fixations allowed (default: 100)

## Output

The `rMASC()` function returns a list containing:

- `results`: Data frame with trial-by-trial results
- `weights`: Vector of attribute weights used
- `parameters`: List of model parameters used
- `raw`: List containing detailed raw data for each trial

## Reference

Gluth, S., Deakin, J., & Rieskamp, J. (2024). A Theory of
Multi-Attribute Search and Choice.
<https://doi.org/10.31234/osf.io/3qzak>

## License

This package is released under the MIT License.
