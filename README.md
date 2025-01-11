
<!-- README.md is generated from README.Rmd. Please edit that file -->

# masc

<!-- badges: start -->
<!-- badges: end -->

## Overview

The `masc` package implements the Multi-Attribute Search and Choice
(MASC) model, a cognitive model of multi-attribute decision making.
Based on the work by Gluth, Deakin, & Rieskamp (2024), this package
provides tools for simulating how people make decisions when comparing
options with multiple attributes.

## Installation

You can install the development version of `masc` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kiantefernandez/masc")
```

## Usage

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

## Contributing

Contributions are welcome! Please file issues or submit pull requests on
the [GitHub repository](https://github.com/kiantefernandez/masc).

## License

This package is released under the MIT License.
