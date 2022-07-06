utils::globalVariables(
  c("tmp_id", "imputation", "imputed_data", "limma_results")
)
#' Run multiple imputation and limma
#'
#' This function is an efficient wrapper that fits the needed gamma regressions,
#' performs multiple imputation and testing with limma (see \code{\link[limma]{limmaUsersGuide}}).
#' It is an efficient wrapper that generates need inputs for imputation and
#' running \code{\link[pair]{run_limma_and_lfc}} with the possibility of using
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
#' @param plot_trend Should the mean-variance trend with the gamma regression be
#'  plotted?
#' @param .robust Input to limma robust option
#'
#' @return A tibble with each imputation as a row. The first column contains the
#'     imputation number, the second contains the imputed data, and the last
#'     column contains the results produced by \code{\link[mavis]{run_limma_and_lfc}}.
#' @export
#'
#' @importFrom utils globalVariables
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_*
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Generate the contrast matrix
#' contrast <- limma::makeContrasts(
#'   contrasts = "ng100-ng50",
#'   levels = design
#' )
#'
#' # Normalize and log-transform the data
#' yeast <- psrn(yeast_prog, "identifier")
#' \dontrun{
#' results <- multiple_imputation_and_limma(yeast, design, contrast, 1000, 5, "identifier", TRUE)
#' }
multiple_imputation_and_limma <- function(data,
                         design,
                         contrast_matrix,
                         imputations,
                         workers = 1,
                         id_col = "id",
                         .robust = TRUE,
                         plot_trend = FALSE) {
  # Fit gamma models
  data <- data %>%
    calculate_mean_sd_trends(design)
  gamma_reg_models <- fit_gamma_regression(data, design)
  gamma_reg_weights <- gamma_reg_models$weights
  if (plot_trend) {
    plot_gamma_regression(data, design, id_col = id_col)
  }
  # Generate imputation input
  LOQ <- data %>%
    purrr::keep(is.numeric) %>%
    unlist(T, F) %>%
    {stats::quantile(., .25, na.rm = T) - 1.5*stats::IQR(., na.rm = T)} %>%
    unname()
  col_order <- names(data)
  missing_data <- data %>%
    dplyr::filter(dplyr::if_any(where(is.numeric), is.na))
  char_cols <- missing_data %>%
    purrr::keep(is.character)
  conditions <- design %>%
    get_conditions()
  impute_nested <- data %>%
    prep_data_for_imputation(conditions, gamma_reg_models$imputation, LOQ)
  # Generate results
  ## Non-missing data
  non_miss_result <- data %>%
    tidyr::drop_na(where(is.numeric)) %>%
    run_limma_and_lfc(design, contrast_matrix, gamma_reg_weights, id_col, .robust = .robust)
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
        "impute_row",
        "char_cols",
        "col_order",
        "run_limma_and_lfc",
        "design",
        "contrast_matrix",
        "gamma_reg_weights",
        "calc_weights",
        "impute_nested",
        "non_miss_result",
        "id_col"
      )
    )
    on.exit(rm(cluster))
  }
  results <- tibble::tibble(
    imputation = seq_len(imputations)
  )
  if (workers != 1) {
    results <- results %>%
      multidplyr::partition(cluster)
  }
  results <- results %>%
    dplyr::mutate(
      # Run imputation
      imputed_data = purrr::map(
        imputation,
        ~ impute(impute_nested, char_cols, col_order)
      ),
      # Run limma
      limma_results = purrr::map(
        imputed_data,
        run_limma_and_lfc,
        design, contrast_matrix, gamma_reg_weights, id_col, NULL, .robust = .robust
      ),
      # Bind non-missing data
      limma_results = purrr::map(
        limma_results,
        dplyr::bind_rows,
        non_miss_result
      )
    )
  if (workers != 1) {
    results <- results %>%
      dplyr::collect() %>%
      dplyr::arrange(imputation)
  }
  return(results)
}
