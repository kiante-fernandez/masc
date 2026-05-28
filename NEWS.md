# masc (development version)

# masc 0.1.0

First public release.

* `rMASC()` simulates the Multi-Attribute Search and Choice (MASC) model of
  Gluth, Deakin and Rieskamp (2026) — sequential information search, Bayesian
  belief updating, and choice.
* New `Sigma_belief` argument enables the multivariate **MASC-C** belief update:
  when the assumed correlation structure is non-diagonal, observing one
  attribute updates beliefs about correlated attributes via a Kalman filter
  ("belief spread"). With a diagonal or `NULL` `Sigma_belief` the model reduces
  exactly to the original univariate MASC update.
* New `Sigma_true` argument generates stimuli with a specified correlation
  structure (a matrix, or a single uniform off-diagonal correlation).
* Bundled `hotelgluth2024` dataset from the hotel-choice experiment.
