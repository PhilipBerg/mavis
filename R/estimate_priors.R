utils::globalVariables(
  c("formula")
)
#' Estimate Gamma priors for sigma
#'
#' @description Estimates the priors for the Bayesian model.
#' @param data A `tibble` or `data.frame` to add gamma priors to
#' @param design_matrix A design matrix for the data (see example)
#' @param gamma_reg A Gamma regression object
#'
#' @return A `tibble` or `data.frame` with the alpha,beta priors
#' @export
#'
#' @examples
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize data
#' yeast <- yeast %>%
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design)
#'
#' # Fit gamma regression
#' gamma_model <- fit_gamma_regression(yeast, sd ~ mean)
#'
#' # Estimate priors
#' gamma_model %>%
#'     estimate_gamma_hyperparameters(yeast, design)
estimate_gamma_hyperparameters <- function(gamma_reg, data, design_matrix){
  UseMethod("estimate_gamma_hyperparameters")
}

#' @export
estimate_gamma_hyperparameters.glm <- function(gamma_reg, data, design_matrix){
  if ("c" %in% colnames(data)) {
    c_col <- rlang::sym("c")
  } else {
    c_col <- NULL
  }
  data %>%
    dplyr::mutate(
      alpha = (1/summary(gamma_reg)$dispersion),
      beta = estimate_beta(gamma_reg, mean, !!c_col, alpha)
    )
}

#' @export
estimate_beta.glm <- function(reg, mean, c, alpha, ...){
  reg_vars <- reg %>%
    formula() %>%
    all.vars()
  reg_vars <- reg_vars[!reg_vars %in% c("mean", "sd")]
  if (length(reg_vars) == 0) {
    reg_vars <- "newdata = data.frame(mean = mean)"
  } else {
    reg_vars <- reg_vars %>%
      paste0(., " = ", ., collapse = ", ") %>%
      paste0("newdata = data.frame(mean = mean, ", ., ")")
  }

  nd <- rlang::parse_expr(reg_vars)
  alpha / rlang::eval_tidy(
    rlang::call2(
      stats::predict.glm,
      object = reg,
      newdata = nd,
      type = 'response'
    )
  )
}
