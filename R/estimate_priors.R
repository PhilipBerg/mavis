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
#' yeast <- yeast_prog %>%
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design)
#'
#' # Fit gamma regression
#' gamma_model <- fit_gamma_regression(yeast, sd ~ mean)
#'
#' # Estimate priors
#' yeast %>%
#'     estimate_gamma_priors(design, gamma_model)
estimate_gamma_priors <- function(data, design_matrix, gamma_reg){
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
