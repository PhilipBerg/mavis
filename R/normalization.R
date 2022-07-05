utils::globalVariables(c("where", "value", "ref", "all_of"))
#' Normalize data to a pseudo-reference
#'
#' This function generates a pseudo-reference by taking the geometric mean of
#' each feature across all samples. Each feature in each sample is then divided
#' by the pseudo-reference. Then, the median ratio of all ratios is used as an
#' estimate to use for normalizing differences in loading concentration. All
#' features in each sample is then divided by their corresponding estimate.
#' All estimates are based on features without missing values.
#' For details see \insertCite{anders2010differential;textual}{mavis}.
#'
#' @param data data.frame
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#' @param load_info logical
#' @param log boolean variable indicating if the data should be log transformed
#'     after normalization
#' @param target target columns to normalize, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns will be used in the normalization if not specified.
#'
#' @return data frame with normalized values if `load_info=FALSE`, if it is `TRUE`
#'    then it returns a list with two tibbles. One tibble containing the
#'    normalized data and one containing the loading info as well as the
#'    estimated normalization factors.
#' @export
#' @importFrom dplyr .data
#' @import utils
#'
#' @examples
#' yeast_psrn <- psrn(yeast_prog, "identifier")
#' @source \url{https://www.nature.com/articles/npre.2010.4282.1}
#' @references
#' \insertAllCited{}
psrn <- function(data,
                 id_col = "id",
                 log = TRUE,
                 load_info = FALSE,
                 target = NULL) {
  target <- rlang::enquo(target)
  target <- check_target(target)
  data_filtered <- data %>%
    tidyr::drop_na(!!target)
  loading_sizes <- calc_loading_size(data_filtered, target)
  pseudo_reference <- data_filtered %>%
    tidyr::pivot_longer(!!target) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      ref = prod(value)^(1 / dplyr::n())
    )
  scaling_factors <- data_filtered %>%
    tidyr::pivot_longer(!!target, names_to = "sample") %>%
    dplyr::left_join(pseudo_reference, by = id_col) %>%
    dplyr::mutate(
      value = value / ref
    ) %>%
    dplyr::select(-ref) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(rle_factor = stats::median(value)) %>%
    dplyr::left_join(loading_sizes, by = "sample")
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <-
      data[scaling_factors$sample[i]] / scaling_factors$rle_factor[i]
  }
  if (log) {
    data <- data %>%
      dplyr::mutate(
        dplyr::across(!!target, log2)
      )
  }
  if (load_info) {
    return(
      list(
        data = data,
        scaling_factors = scaling_factors
      )
    )
  } else {
    return(data)
  }
}

utils::globalVariables(c("load_size", "lfc", "A", "w"))
#' Normalization by Trimmed m Means
#'
#' Method for trimmed m means normalization. It is based on the method described
#' in \insertCite{robinson2010scaling;textual}{mavis}. Though, instead of library
#' size as was used in the original method, here we use the loading size which
#' we define as the sum of all features. If no reference sample is specified, it
#' uses the sample with the lowest coefficient of variation as default.
#' All estimates are based on features without missing values.
#'
#' @param data data.frame containing the data to normalize
#' @param trim_M percent of fold-change values to trim
#' @param trim_A percent of means to trim
#' @param log Return log2 transformed values?
#' @param load_info Return loading info?
#' @param target target columns to normalize, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns will be used in the normalization if not specified.
#' @param reference_sample  Specify a reference sample to normalize to, if not
#'     provided, the sample with the lowest coefficient of variation will be used
#'
#' @return data frame with normalized values if `load_info=FALSE`, if it is `TRUE`
#'    then it returns a list with two tibbles. One tibble containing the
#'    normalized data and one containing the loading info as well as the
#'    estimated normalization factors.
#' @export
#' @import utils
#'
#' @examples
#' tmm(yeast_prog)
#' @source \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25}
#' @references
#' \insertAllCited{}
tmm <- function(data,
                trim_M = .3,
                trim_A = .05,
                log = TRUE,
                load_info = FALSE,
                target = NULL,
                reference_sample = NULL) {
  target <- rlang::enquo(target)
  target <- check_target(target)
  data_filtered <- data %>%
    tidyr::drop_na(!!target, reference_sample)
  if (is.null(reference_sample)) {
    reference_sample <- calc_cv(data_filtered, target)
    reference_sample <-
      names(reference_sample)[reference_sample == min(reference_sample)]
  }
  loading_sizes <- calc_loading_size(data_filtered, target) %>%
    dplyr::bind_rows(
      calc_loading_size(data_filtered, reference_sample)
    ) %>%
    dplyr::distinct()
  reference_loading_size <-
    loading_sizes$load_size[loading_sizes$sample == reference_sample]
  reference_sample <- rlang::sym(reference_sample)
  scaling_factors <- data_filtered %>%
    dplyr::select(!!target, reference_sample) %>%
    tidyr::pivot_longer(-!!reference_sample, names_to = "sample") %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::mutate(
      w = (load_size - value) / (load_size * value) +
        (reference_loading_size - !!reference_sample) /
          (reference_loading_size * !!reference_sample),
      lfc = log2((value / load_size) /
        (!!reference_sample / reference_loading_size)),
      A = .5 * log2((value / load_size) *
        (!!reference_sample / reference_loading_size))
    ) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(
      dplyr::between(
        lfc, stats::quantile(lfc, trim_M), stats::quantile(lfc, 1 - trim_M)
      ) &
        dplyr::between(
          A, stats::quantile(A, trim_A), stats::quantile(A, 1 - trim_A)
        )
    ) %>%
    dplyr::summarise(
      tmm_factor = 2^(sum(lfc * w) / sum(w))
    ) %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::add_row(
      sample = as.character(reference_sample),
      tmm_factor = 1,
      load_size = reference_loading_size
    )
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <- data[scaling_factors$sample[i]] /
      scaling_factors$tmm_factor[i]
  }
  if (log) {
    data <- data %>%
      dplyr::mutate(
        dplyr::across(!!target, log2)
      )
  }
  if (load_info) {
    return(
      list(
        data = data,
        scaling_factors = scaling_factors
      )
    )
  } else {
    return(data)
  }
}

calc_cv <- function(data, targets) {
  data %>%
    dplyr::select(!!targets) %>%
    purrr::map_dbl(
      ~ sd(.x) / mean(.x)
    )
}

calc_loading_size <- function(data, targets) {
  data %>%
    dplyr::select(!!targets) %>%
    colSums() %>%
    tibble::enframe(name = "sample", value = "load_size")
}
