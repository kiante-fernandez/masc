#' masc: Multi-Attribute Search and Choice Model
#'
#' @description
#' Implementation of the Multi-Attribute Search and Choice (MASC) model for
#' decision making. This model simulates how people make decisions when comparing
#' options with multiple attributes by sequentially sampling information about
#' different attributes until reaching a decision.
#'
#' @details
#' The main function provided by this package is:
#' \itemize{
#'   \item \code{\link{rMASC}}: Generate choices using the MASC model
#' }
#'
#' @references
#' Gluth, S., Deakin, J., & Rieskamp, J. (2024). A Theory of Multi-Attribute Search
#' and Choice.
#'
#' @docType package
#' @name masc-package
#' @keywords internal
"_PACKAGE"

#' @useDynLib masc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats pnorm qnorm rbeta rmultinom rnorm
NULL
