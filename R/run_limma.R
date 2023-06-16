utils::globalVariables(c("logFC", "AveExpr", "P.Value", "rowid"))
#' Run limma and calculate LFC
#'
#' Function to calculate LFC, run limma using \code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}},
#' and \code{\link[limma]{eBayes}} with the flag robust = TRUE.
#' If neither `weights` or `gamma_reg_model` is provided, then
#' \code{\link[limma]{lmFit}} will be run without precision weights and \code{trend = TRUE} (i.e., limma-trend).
#' If both are provided, then it will default to using the `gamma_reg_model`
#' model to produce the weights.
#'
#' @param data a `data.frame` with the samples and feature ids.
#' @param design a design or model matrix produced by
#'     \code{\link[stats]{model.matrix}}.
#' @param contrast_matrix a contrast matrix produced by
#'     \code{\link[limma]{makeContrasts}}.
#' @param gamma_reg_model the regression model produced by
#'     `fit_gamma_regression` (see \code{\link[mavis]{fit_gamma_regression}}) or any
#'      `glm` with `formula` `sd ~ ...`.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'
#' @param weights a matrix of precision weights.
#' @param .robust Logical value if limma::eBayes should use robust estimator (TRUE) or not (FALSE)
#'
#' @return a `tibble` with the id_col, then one p_val_* and lfc_* for each
#'     comparison (*) in the contrast matrix
#' @export
#'
#' @importFrom rlang :=
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
#' yeast <- psrn(yeast, "identifier")
#'
#' # Fit the gamma regressions
#' gamma_reg_model <- yeast %>%
#'   calculate_mean_sd_trends(design) %>%
#'   fit_gamma_regression(formula = sd ~ mean)
#'
#' # Exemplify on the non-missing data
#' yeast <- tidyr::drop_na(yeast)
#'
#' results <- run_limma(
#'   yeast,
#'   design,
#'   contrast,
#'   gamma_reg_model,
#'   "identifier"
#' )
run_limma <- function(data,
                              design,
                              contrast_matrix,
                              gamma_reg_model = NULL,
                              id_col = "id",
                              weights = NULL,
                              .robust = TRUE) {
  row_names <- data[[id_col]]
  condi <- design %>%
    get_conditions()
  data <- data %>%
    dplyr::select(tidyr::matches(condi), dplyr::any_of('c')) %>%
    as.data.frame()
  rownames(data) <- row_names
  # Run LIMMA
  if (!is.null(gamma_reg_model)) {
    weights <- calc_weights(data, id_col, design, gamma_reg_model)
    rownames(weights) <- row_names
    trend <- FALSE
  } else if (!is.null(weights)) {
    if((dim(weights) != dim(data))){
      msg <- "Incorect dimensions between data and weights."
      incorrect_dims <- which(dim(weights) != dim(data))
      for (i in seq_along(incorrect_dims)) {
        msg[i + 1] <- dplyr::case_when(
          incorrect_dims[i] == 1 ~ glue::glue("Data has {dim(data)[incorrect_dims[i]]} rows while weights have {dim(weights)[incorrect_dims[i]]}"),
          incorrect_dims[i] == 2 ~ glue::glue("Data has {dim(data)[incorrect_dims[i]]} columns while weights have {dim(weights)[incorrect_dims[i]]}")
        )
      }
      rlang::abort(
        stringr::str_flatten(msg, "\n")
      )
    }
    trend <- FALSE
  } else{
    trend <- TRUE
  }
  data <- data %>%
    dplyr::select(dplyr::matches(condi)) %>%
    as.matrix()
  hits <- limma::lmFit(data, design, weights = weights) %>%
    limma::contrasts.fit(contrast_matrix) %>%
    limma::eBayes(robust = .robust, trend = trend)
  # Extract p-value and LFC from comparisons
  colnames(contrast_matrix) %>%
   stats:: setNames(colnames(contrast_matrix)) %>%
    purrr::map(
      ~limma::topTable(hits, .x, number = Inf, adjust.method = 'none') %>%
        tibble::rownames_to_column(id_col) %>%
        tibble::as_tibble()
    ) %>%
    purrr::imap(
      dplyr::mutate
    ) %>%
    purrr::map(
      dplyr::rename, comparison = dplyr::last_col()
    ) %>%
    dplyr::bind_rows() %>%
    dplyr::select(
      all_of(id_col),
      lfc = logFC,
      mean = AveExpr,
      p_val = P.Value,
      comparison
    ) %>%
    dplyr::mutate(
      comparison = stringr::str_replace(comparison, '-', ' vs ')
    )
}
