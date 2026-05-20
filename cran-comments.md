## Submission

This is a new submission: the first release of `masc`, version 0.1.0.

## R CMD check results

Local `R CMD check --as-cran` (macOS, R 4.x) and win-builder (R-devel and
R-release):

    0 errors | 0 warnings | 0 notes

On CRAN's incoming checks a "New submission" NOTE is expected for a first
release.

The DESCRIPTION may trigger a "Possibly mis-spelled words" NOTE for the
following, which are all intentional (a model name, the authors of the cited
work, and standard technical terms): MASC, multiattribute, Bayesian, Kalman,
Gluth, Deakin, Rieskamp.

The Description field cites the implemented model via DOI
<doi:10.1037/rev0000614>, which resolves to the published article.

## Test environments

* local macOS, R 4.x (R CMD check --as-cran)
* win-builder, R-devel and R-release (devtools::check_win_devel / _release)

## revdepcheck results

There are no downstream dependencies (new package).
