#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point using a Gamma regression.
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
  if(!is.null(data$c)){
    pred <- ~stats::predict.glm(
      gam_reg,
      newdata = data.frame(mean = ., c = c),
      type = "response"
    )
  } else {
    pred <- ~predict.glm(
      gam_reg,
      newdata = data.frame(mean = .),
      type = "response"
    )
  }
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')
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
