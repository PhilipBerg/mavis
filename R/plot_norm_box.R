utils::globalVariables(c("where", "value", "name", "method", "med", "condition"))
#' Generates a box-plot for `log2` transformed raw, tmm-, and psn- normalized
#'     values.
#'
#' This function can be used to produce a visual aid to for selection
#' normalization method. Currently \code{\link[mavis]{tmm}} and \code{\link[mavis]{psrn}}
#' are used. Ideally, after the data is normalized all the boxplots should have
#' their median aligned close to the global trend line.
#'
#' @param data data.frame containing the data to normalize
#' @param id_col a character for the name of the column containing the
#'      name of the features in data (e.g., peptides, proteins, etc.)
#' @param trim_M percent of fold-change values to trim
#' @param trim_A percent of means to trim
#' @param reference_sample Specify a reference sample to normalize to in the
#'     \code{\link[mavis]{tmm}} method
#' @param norm_target target columns to normalize, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns will be used in the normalization if not specified.
#' @param plot_target target columns to plot, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns are plotted.
#'
#' @return ggplot with samples on the x-axis and observed values on the y-axis,
#'     different colors correspond to the raw or normalized data. Lines
#'     corresponds to the global median across all samples.
#' @export
#'
#' @examples
#' plot_norm_box(yeast, "identifier")
#' # Plot only ng50 samples
#' plot_norm_box(yeast, "identifier", plot_target = contains("ng50"))
plot_norm_box <- function(data,
                          id_col = "id",
                          trim_M = .3,
                          trim_A = .05,
                          norm_target = NULL,
                          plot_target = NULL,
                          reference_sample = NULL) {
  norm_target <- rlang::enquo(norm_target)
  norm_target <- check_target(norm_target)
  plot_target <- rlang::enquo(plot_target)
  plot_target <- check_target(plot_target)
  psrn <- data %>%
    baldur::psrn(
      id_col = id_col,
      target = !!norm_target
    ) %>%
    dplyr::select(id_col, !!plot_target) %>%
    dplyr::rename_with(~ paste0(., "_psrn"), where(is.numeric))
  tmm <- data %>%
    tmm(
      trim_M = trim_M,
      trim_A = trim_A,
      target = !!norm_target,
      reference_sample = reference_sample
    ) %>%
    dplyr::select(id_col, !!plot_target) %>%
    dplyr::rename_with(~ paste0(., "_tmm"), where(is.numeric))
  data %>%
    dplyr::select(id_col, !!plot_target) %>%
    dplyr::mutate(
      dplyr::across(!!plot_target, log2)
    ) %>%
    dplyr::rename_with(~ paste0(., "_raw"), where(is.numeric)) %>%
    dplyr::left_join(psrn, by = id_col) %>%
    dplyr::left_join(tmm, by = id_col) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    tidyr::extract(name, c("condition", "method"), "^(.*)_(.*)$") %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(
      med = stats::median(value, na.rm = T)
    ) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(condition, value, fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = med, color = method)) +
    ggplot2::theme_bw()
}
