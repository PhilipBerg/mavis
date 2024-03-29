utils::globalVariables(
  c("tmp_id", "imputation", "imputed_data", "limma_results")
)
#' Run multiple imputation and limma
#'
#' This function is an efficient wrapper that fits the needed gamma regressions,
#' performs multiple imputation and testing with limma (see \code{\link[limma]{limmaUsersGuide}}).
#' It is an efficient wrapper that generates need inputs for imputation and
#' running \code{\link[mavis]{run_limma}} with the possibility of using
#' \code{\link[multidplyr]{multidplyr-package}} to paralellize the computation.
#' It also calls produces the mean-variance trends if `plot_trend` is `TRUE`.
#'
#'
#' @param data The data to run the pipeline on, missing values should have NA
#'     values.
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param contrast_matrix A contrast matrix of comparisons to perfrom see
#'     \code{\link[limma]{makeContrasts}} for details.
#' @param imputations Number of imputations to perfrome.
#' @param workers Number of workers (processes) to run the pipeline with.
#'    Any value >1 will run the pipeline with parallel computing using the
#'    \code{\link[multidplyr]{multidplyr-package}}.
#' @param id_col A character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#' @param .robust Input to limma robust option
#' @param weights logical value if weights should be used.
#' @param plot_trend Should the mean-variance trend with the gamma regression be
#'  plotted?
#' @param formula_imputation Formula for the regression model used for imputation
#' @param formula_weights Formula for the regression model used for weights
#' @param imp_pars Imputation parameters to use; as inferred by [get_imputation_pars].
#' @param ... Additional arguments to \code{\link[mavis]{multi_trend_partitioning}}
#'
#' @return A tibble with each imputation as a row. The first column contains the
#'     imputation number, the second contains the imputed data, and the last
#'     column contains the results produced by \code{\link[mavis]{run_limma}}.
#' @export
#'
#' @importFrom utils globalVariables
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important.
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Generate the contrast matrix
#' contrast <- limma::makeContrasts(
#'   contrasts = "ng100-ng50",
#'   levels = design
#' )
#'
#' # Normalize and log-transform the data
#' yeast <- psrn(yeast, "identifier")
#' \donttest{
#' results <- multiple_imputation_and_limma(
#'   yeast,
#'   design,
#'   contrast,
#'   10, # Just running 10 imputation for the example purpose, but can be large > 1000
#'   1,  # Number of parallel workers to use (will drastically increase the speed)
#'   "identifier"
#' )
#' }
multiple_imputation_and_limma <- function(data,
                                          design,
                                          contrast_matrix,
                                          imputations,
                                          workers = 1,
                                          id_col = "id",
                                          .robust = TRUE,
                                          weights = TRUE,
                                          plot_trend = FALSE,
                                          formula_imputation = sd_p ~ mean,
                                          formula_weights = sd ~ mean,
                                          imp_pars = NULL,
                                          ...) {
  # Fit gamma models
  data <- data %>%
    calculate_mean_sd_trends(design)
  if (
    !'c' %in% names(data) &
    any(stringr::str_detect(
        as.character(formula_weights),
         'c')
    )
  ) {
    data <- data %>%
      multi_trend_partitioning(design, formula_weights, ...)
  }
  if (weights) {
    gamma_reg_weights <- fit_gamma_regression(data, formula_weights)
  } else{
    gamma_reg_weights <- NULL
  }
  if (plot_trend) {
    plot_gamma_regression(data, design, ...)
  }
  # Generate imputation input
  if (
    !'c' %in% names(data) &
    any(stringr::str_detect(
      as.character(formula_imputation),
      'c')
    ) & is.null(imp_pars)
  ) {
    cat('Trend paritioning for imputation\n')
    imp_data <- data %>%
      multi_trend_partitioning(design, formula_imputation, ...)
  } else {
    imp_data <- data
  }

  if (is.null(imp_pars)) {
    imp_pars <- imp_data %>%
      get_imputation_pars(design, formula_imputation, workers, ...)
  }

  if (workers != 1) {
    cluster <- multidplyr::new_cluster(workers)
    multidplyr::cluster_library(
      cluster,
      c(
        "dplyr",
        "stringr",
        "tidyr",
        "purrr",
        "limma",
        "tibble"
      )
    )
    multidplyr::cluster_copy(
      cluster,
      c(
        "impute",
        "run_limma",
        "design",
        "contrast_matrix",
        "gamma_reg_weights",
        "calc_weights",
        "imp_pars",
        "id_col",
        "data"
      )
    )
    on.exit(rm(cluster))
  }
  results <- tibble::tibble(
    imputation = seq_len(imputations)
  ) %>%
    dplyr::mutate(
      # Run imputation
      imputed_data = purrr::map(
        imputation,
        ~ impute(data, imp_pars)
      )
    )
  if (workers != 1) {
    results <- results %>%
      multidplyr::partition(cluster)
  }
  results <- results %>%
    dplyr::mutate(
      # Run limma
      limma_results = purrr::map(
        imputed_data,
        run_limma,
        design, contrast_matrix, gamma_reg_weights, id_col, NULL, .robust = .robust
      )
    )
  if (workers != 1) {
    results <- results %>%
      dplyr::collect() %>%
      dplyr::arrange(imputation)
  }
  return(results)
}
