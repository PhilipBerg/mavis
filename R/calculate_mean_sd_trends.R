#' Calculate the Mean-Variance trend
#'
#' @description Calculates the mean and sd of the rows.
#'
#' @param data A `tibble` or `data.frame` to annotate with mean and sd
#' @param design_matrix A design matrix for the data (see example)
#' @param sdev Which centered second order moment to calculate. Can be "sample"
#'   for the sample standard deviation, "pooled" for the pooled standard
#'   deviation, or can be both of them to calculate both.
#'
#' @return A `tibble` or `data.frame` with the mean and sd vectors
#' @export
#'
#' @examples
#' #Example
calculate_mean_sd_trends <- function(data, design_matrix, sdev = c("sample", "pooled")){
  sdev <- match.arg(sdev, c("sample", "pooled"), several.ok = TRUE)
  conditions <- design_matrix %>%
    colnames() %>%
    paste0(collapse = '|')
  data <- data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions)), na.rm = T)
    )
  if (length(sdev) == 2) {
    data %>%
      dplyr::mutate(
        sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd, na.rm = T)
      ) %>%
      pooled_sd(design_matrix)
  } else if (sdev == 'sample') {
    data %>%
      dplyr::mutate(
        sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd, na.rm = T)
      )
  } else {
    data %>%
      pooled_sd(design_matrix)
  }
}

utils::globalVariables(
  c("sd_p")
)
pooled_sd <- function(data, design) {
  condi  <- get_conditions(design)

  weights <- data %>%
    dplyr::select(matches(condi)) %>%
    split.default(stringr::str_extract(colnames(.), condi)) %>%
    purrr::map(
      ~ apply(.x, 1, \(x) sum(is.na(x)))
    ) %>%
    purrr::map2(colSums(design),
                ~ .y - .x - 1
    ) %>%
    purrr::map(
      ~ dplyr::if_else(.x>0, .x, 0)
    ) %>%
    tibble::as_tibble()

  data %>%
    dplyr::select(matches(condi)) %>%
    split.default(stringr::str_extract(colnames(.), condi)) %>%
    purrr::map2(weights,
                ~ apply(.x, 1, var, na.rm = T)*.y
    ) %>%
    tibble::as_tibble() %>%
    dplyr::reframe(
      sd_p = apply(dplyr::across(everything()), 1, sum, na.rm = T)/rowSums(weights),
      sd_p = sqrt(sd_p)
    ) %>%
    dplyr::bind_cols(
      dplyr::select(data, -dplyr::any_of('sd_p')), .
    )
}
