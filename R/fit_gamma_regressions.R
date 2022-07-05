utils::globalVariables(c(".", "sd", "model"))
#' Functions for Fitting the Mean-Variance Gamma Regression Models
#'
#' `fit_gamma_regressions` is a wrapper function that calls both
#'    \code{fit_gamma_imputation} and
#'    \code{fit_gamma_weights}. It returns a list containing
#'    the models for imputation in `$imputation` and the weights in `$weights`.
#'    `fit_gamma_imputation` returns a list named according to the different
#'    conditions and `fit_gamma_weights` returns a `glm` object containing the
#'    gamma regression for the mean-variance trend.
#'
#' @import utils
#'
#' @param data a `data.frame` to generate the mean-variance trends for. It
#'     should contain columns with conditions named as the column names in
#'     `design` (presumably with some suffix).
#' @param design a design or model matrix as produced by
#'  \code{\link[stats]{model.matrix}} with column names corresponding to the
#'  different conditions.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'
#' @return `fit_gamma_imputation` returns a named
#'     list where the names corresponds to the conditions. Each index contains
#'     a `glm` object with the gamma regression for the mean-variance trend.
#'     `fit_gamma_weights` returns a `glm` object with the gamma regression
#'     for the precision weights. `fit_gamma_regressions` returns a list with
#'     the results from `fit_gamma_imputation` in `$imputation` and the results
#'     from `fit_gamma_weights` in `$weights`.
#' @name Mean-Variance_Gamma_Regressions
NULL

#' @describeIn Mean-Variance_Gamma_Regressions Wrapper function that runs both
#'     `fit_gamma_imputation` and `fit_gamma_weights`
#' @export
#'
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_*
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log-transform the data
#' yeast <- psrn(yeast_prog, "identifier")
#'
#' # Fit all gamma regression models for the mean-variance trends
#' all_gamma_models <- fit_gamma_regressions(yeast, design, "identifier")
fit_gamma_regressions <- function(data, formula = sd ~ mean + c) {
    glm(formula, Gamma(log), data)
}

#' @describeIn Mean-Variance_Gamma_Regressions Function that generates per
#'     condition mean-variance trends used in the imputation procedure.
#'     Each id in the `id_col` gets one mean and variance calculated for each
#'     condition. One model is then fitted per condition.
#' @return gamma regressions for each condition
#' @export
#'
#' @examples
#'
#' # Fit the gamma regression models for the mean-variance trend used in the
#' # imputation procedure
#' gamma_imputation_models <- fit_gamma_imputation(yeast, design, "identifier")
fit_gamma_imputation <- function(data, design, id_col = "id") {
  data %>%
    prep_data_for_gamma_imputation_regression(design, id_col) %>%
    dplyr::group_by(name) %>%
    fit_gamma_model() %>%
    split.data.frame(.$name) %>%
    extract_model()
}

#' @describeIn Mean-Variance_Gamma_Regressions Function to produce the
#'     mean-variance trend used to calculate the precision weights used in
#'     \code{\link[limma]{lmFit}}. Each id in the `id_col` gets one mean and one
#'     variance across all conditions and one model is then fitted for all
#'     mean-variance pairs.
#'
#' @return gamma regression for weights
#' @export
#'
#' @examples
#'
#' # Fit the gamma regression model for the mean-variance trend used for
#' # estimating the precision weights used in limma
#' gamma_weight_model <- fit_gamma_weights(yeast, design, "identifier")
fit_gamma_weights <- function(data, design, id_col = "id") {
  data %>%
    prep_data_for_gamma_weight_regression(design, id_col) %>%
    fit_gamma_model() %>%
    extract_model()
}

prep_data_for_gamma_weight_regression <- function(data, design, id_col = "id") {
  data %>%
    pivot_data_for_gamma_regression(design) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    calc_mean_sd_trend()
}

calc_mean_sd_trend <- function(data) {
  data %>%
    dplyr::summarise(
      mean = mean(value, na.rm = T),
      sd = stats::sd(value, na.rm = T),
      .groups = "drop"
    ) %>%
    dplyr::filter(sd > 0)
}

fit_gamma_model <- function(data) {
  data %>%
    dplyr::summarise(
      model = list(stats::glm(sd ~ mean, stats::Gamma(log))),
      .groups = "drop"
    )
}

prep_data_for_gamma_imputation_regression <- function(data,
                                                      design,
                                                      id_col = "id") {
  conditions <- design %>%
    get_conditions()
  data %>%
    pivot_data_for_gamma_regression(design) %>%
    dplyr::mutate(
      name = stringr::str_replace(
        name,
        paste0("(", conditions, ")", ".*"), "\\1"
      )
    ) %>%
    dplyr::group_by(name, .data[[id_col]]) %>%
    dplyr::filter(sum(!is.na(value)) >= 2) %>%
    calc_mean_sd_trend()
}

pivot_data_for_gamma_regression <- function(data, design) {
  conditions <- design %>%
    get_conditions()
  data %>%
    tidyr::pivot_longer(dplyr::matches(conditions))
}

extract_model <- function(data) {
  if (is.data.frame(data)) {
    data %>%
      magrittr::use_series(model) %>%
      magrittr::extract2(1)
  } else if (is.list(data)) {
    data %>%
      purrr::map(magrittr::use_series, model) %>%
      purrr::map(magrittr::extract2, 1)
  } else {
    error_message <- paste(
      "Data of class",
      paste0(class(data), collapse = ", "),
      "not suppported."
    )
    rlang::abort(error_message)
  }
}
