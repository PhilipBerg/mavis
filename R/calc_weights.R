#' Calculate Precision Weights From Gamma Regression
#'
#' Calculates the precision weights for each measurement.
#'
#' \code{calc_weights} takes as input a data frame and a \code{\link[stats]{glm}}  object produced
#' by fit_gamma_regression,
#' see \code{\link[mavis]{fit_gamma_regression}} for details.
#' For all numeric columns, it predicts the standard deviation using the gamma
#' regression. It then squares and takes the reciprocal of each value to generate the
#' precision weights.
#'
#' @param data A `data.frame` who's quantitative columns should be converted to
#'  precision weights.
#' @param gamma_reg_model a `glm` object as produced by
#'   see \code{\link[mavis]{fit_gamma_regression}} for details.
#'
#' @param identifier A string for the unique identifiers in the data
#' @param design A design matrix for the data
#'
#' @return The same `data.frame` but with all quantitative values replaced by
#'  their precision weights.
#'
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_weights
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log-transform the data
#' yeast <- psrn(yeast, "identifier") %>%
#'   # Add row means and variances
#'   calculate_mean_sd_trends(design)
#'
#' # Fit the gamma regression model for the mean-variance trend
#' gamma_model <- fit_gamma_regression(yeast, sd ~ mean)
#'
#' # Generate the weights for the yeast data
#' calc_weights(yeast, 'identifier', design, gamma_model)
calc_weights <- function(data, identifier, design, gamma_reg_model) {
  estimate_uncertainty(data, identifier, design, gamma_reg_model) %>%
    magrittr::raise_to_power(-2)
}
