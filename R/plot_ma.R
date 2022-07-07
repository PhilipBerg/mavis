utils::globalVariables(c("significant", "median_lfc"))
#' Generate a MA-plot of the analysis.
#'
#' @param hits a `tibble` produced by \code{\link[mavis]{extract_results}}
#' @param data the `data.frame` used in the analysis.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[mavis]{extract_results}}.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[mavis]{extract_results}}.
#' @param alpha The alpha cut-off for considering a p-value significant.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[mavis]{extract_results}}.
#' @param abs_lfc If a LFC threshold should also be used in the decision.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[mavis]{extract_results}}.
#'
#' @return a `ggplot2` of the distribution of the hits.
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
#'
#' results <- run_pipeline(yeast, design, contrast, 1000, 5, "identifier", TRUE)
#' imputation_summary <- extract_results(yeast, results, .05, 1, "fdr", "identifier")
#' plot_ma(imputation_summary)
#' }
plot_ma <- function(hits, data = NULL, id_col = NULL, alpha = NULL, abs_lfc = NULL) {
  if (is.null(data)) {
    hits <- hits %>%
      dplyr::mutate(
        significant = dplyr::case_when(
          binom_p_value < .05 & median_lfc < 0 ~ "Down",
          binom_p_value < .05 & median_lfc > 0 ~ "Up",
          T ~ "Not significant"
        ),
        significant = factor(significant,
                             levels = c("Not significant", "Up", "Down")
        )
      )
  }else {
    hits <- hits %>%
      dplyr::mutate(
        significant = dplyr::case_when(
          p_val < alpha & lfc > abs_lfc ~ "Up",
          p_val < alpha & lfc < -abs_lfc ~ "Down",
          T ~ "Not significant"
        ),
        significant = factor(significant,
                             levels = c("Not significant", "Up", "Down")
        )
      ) %>%
      dplyr::rename(median_lfc = lfc, median_mean = mean)
  }
  hits %>%
    dplyr::arrange(significant) %>%
    ggplot2::ggplot(
      ggplot2::aes(median_mean, median_lfc, color = significant)
    ) +
    ggplot2::geom_point(size = 1 / 2) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(comparison ~ .) +
    ggplot2::xlab(expression(bar(Y))) +
    ggplot2::ylab("LFC")
}
