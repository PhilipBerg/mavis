#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point using
#'  a Gamma regression.
#' @param data A `tibble` or `data.frame`
#' @param identifier Id column
#' @param design_matrix Cell mean design matrix for the data
#' @param gam_reg A gamma regression model
#'
#' @return A matrix with the uncertainty
#' @export
#'
#' @examples # Please see vignette('baldur') for examples
estimate_uncertainty <- function(data, identifier, design_matrix, gam_reg){
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')
  pred <- get_pred_fun(gam_reg)
  data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric),
                    !!pred
      )
    ) %>%
    dplyr::select(dplyr::matches(condi_regex)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[identifier]])
}

get_pred_fun <- function(reg) {
  model <- rlang::enquo(reg) %>%
    rlang::as_name() %>%
    rlang::sym()
  reg_vars <- reg %>%
    stats::formula() %>%
    all.vars()
  reg_vars <- reg_vars[!reg_vars %in% c('mean', 'sd', 'sd_p')]
  if (length(reg_vars) != 0) {
    reg_vars <- reg_vars %>%
      paste0(., ' = ', ., collapse = ', ') %>%
      paste0("mean = .x, ", .)
  } else {
    reg_vars <- "mean = .x"
  }
  reg_vars <- paste0("data.frame(",  reg_vars, ")")
  nd <- rlang::parse_expr(reg_vars)
  rlang::quo(
    ~ stats::predict.glm(
      object = !!model,
      newdata = !!nd,
      type = "response"
    )
  )
}
