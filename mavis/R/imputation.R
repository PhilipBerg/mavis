#' Single imputation
#'
#' Performs a single imputation run and returns the data with NA values replaced
#' by imputed values.
#'
#' @param data a `data.frame` to perform the imputation on, missing values should
#' be `NA`.
#' @param design a design or model matrix as produced by
#'  \code{\link[stats]{model.matrix}} with column names corresponding to the
#'  different conditions.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'
#' @return a `data.frame` with `NA` values replaced by imputed values.
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_weights
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_prog %>%
#'   # Normalize and log-transform the data
#'   psrn("identifier") %>%
#'   # Run the imputation
#'   single_imputation(design, "identifier")
single_imputation <- function(data,
                              design,
                              id_col = "id",
                              gam_reg) {
  conditions <- design %>%
    get_conditions()
  LOQ <- data %>%
    dplyr::select(dplyr::matches(conditions)) %>%
    unlist(T, F) %>%
    {stats::quantile(., .25, na.rm = T) - 1.5*stats::IQR(., na.rm = T)} %>%
    unname()
  order <- data %>%
    colnames()
  data <- data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions)), na.rm = T)
    ) %>%
    prep_data_for_imputation(conditions, gam_reg, LOQ)
  aux_cols <- data[[1]] %>%
    dplyr::select(-c(mean_condi, data))
  data %>%
    purrr::map(
      dplyr::select, c(sd, mean_condi, data)
    ) %>%
    impute(order) %>%
    purrr::map(
      dplyr::select, -c(sd, mean_condi)
    ) %>%
    dplyr::bind_cols(aux_cols, .) %>%
    dplyr::select(dplyr::all_of(order))
}


impute <- function(data, order) {
  data %>%
    purrr::map(
      dplyr::mutate,
      data = purrr::pmap(list(mean_condi, sd, data), impute_row)
    ) %>%
    purrr::map(tidyr::unnest, data)
}

impute_nest <- function(data, condition, gamma_reg_model, LOQ) {
  if (anyNA(data)) {
    data <- data %>%
      dplyr::mutate(
        mean_condi = rowMeans(dplyr::across(dplyr::contains(condition)), na.rm = T),
        mean_condi = tidyr::replace_na(mean_condi, LOQ)
      )
    data[["sd"]] <- stats::predict(gamma_reg_model, data, type = "response")
    data %>%
      tidyr::nest(data = dplyr::contains(condition))
  } else {
    return(data)
  }
}

impute_row <- function(mean_condi, sd, data) {
  if (!anyNA(data)) {
    return(data)
  } else {
    data <- as.data.frame(data)
    data[is.na(data)] <- stats::rnorm(n = sum(is.na(data)), mean = mean_condi, sd = sd)
    return(data)
  }
}


prep_data_for_imputation <- function(data, conditions, gamma_reg_imputation, LOQ) {
  tmp_cols <- data %>%
    dplyr::select(-dplyr::matches(conditions))
  data %>%
    split.default(stringr::str_extract(names(.), conditions)) %>%
    purrr::map(dplyr::bind_cols, tmp_cols) %>%
    purrr::imap(impute_nest, gamma_reg_imputation, LOQ)
}

estimate_loq <- function(data) {
  data %>%
    purrr::keep(is.numeric) %>%
    unlist(T, F) %>%
    {
      stats::quantile(., .25, na.rm = T) - 1.5 * stats::IQR(., na.rm = T)
    } %>%
    unname()
}
