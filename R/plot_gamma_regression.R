#'@importFrom Rdpack reprompt
utils::globalVariables(c(".", "sd", "model"))
#' Function for plotting the mean-variance gamma regressions
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the precision weights and imputation.
#'
#' @name plt_gam
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param ... Additional arguments to \code{\link[mavis]{multi_trend_partitioning}}.
#' @param sd_type To plot the sample standard deviation ("sd") or pooled estimator
#' ("sd_p").
#'
#' @return a plot with the mean-variance trend used for the precision
#' weights on the left side, and the trend lines used for the imputation on the
#' right side.
#' @export
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- psrn(yeast, "identifier")
#'
#' # Generate the plots
#' \dontrun{
#' plot_gamma_regression(yeast, design, verbose = FALSE)
#' }
plot_gamma_regression <- function(data, design, sd_type = 'sd', ...) {
  data <- data %>%
    calculate_mean_sd_trends(design)
  base_plot <- data %>%
    plot_gamma(sd_type)
  part_plot <- data %>%
    plot_gamma_partition(design, ...)
  plots <- cowplot::plot_grid(base_plot, part_plot)
  title <- cowplot::ggdraw() +
    cowplot::draw_label("Mean-Variance trends", fontface = "bold")
  cowplot::plot_grid(
    title,
    plots,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
}

plot_mean_sd_trend <- function(data, sd_type) {
  y_tit <- set_y_tit(sd_type)
  data %>%
    ggplot2::ggplot(ggplot2::aes(mean, !!rlang::sym(sd_type))) +
    ggplot2::geom_point(size = 1 / 10) +
    ggplot2::geom_smooth(
      method = stats::glm,
      formula = y ~ x,
      method.args = list(family = stats::Gamma(log)),
      fullrange = TRUE,
      se = F,
      color = 'grey'
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = expression(bold(bar(y))), y = y_tit
    )
}

#' Function for plotting the gamma regression used for mean-variance normalization
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the precision weights
#'
#' @inheritParams plt_gam
#'
#' @return a plot with the mean-variance trend used for the precision
#' weights.
#' @export
#' @importFrom tidyr drop_na
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- psrn(yeast, "identifier")
#'
#' # Generate the plot
#' yeast %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma()
plot_gamma <- function(data, sd_type = 'sd') {
  data %>%
    tidyr::drop_na(sd_type) %>%
    plot_mean_sd_trend(sd_type) +
    ggplot2::ggtitle("Before Partitioning")
}

#' Function for plotting the gamma regression used for imputation
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the imputation.
#' @inheritParams plt_gam
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param score If `NULL` does nothing, if a real value, annotates the plot with
#' the `score`.
#' @param ... Additional arguments to \code{\link[mavis]{trend_partitioning}}.
#'
#' @return a plot with the mean-variance trend used for the trend lines used in
#'  the imputation.
#' @export
#'
#' @importFrom grDevices adjustcolor
#' @importFrom stats setNames
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- psrn(yeast, "identifier")
#'
#' # Generate the plot
#' \dontrun{
#' yeast %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma_partition(design, verbose = FALSE)
#' }
plot_gamma_partition <- function(data, design, score = NULL, ...) {
  if(!'c' %in% names(data)){
    data <- data %>%
      multi_trend_partitioning(design, ...)
  }
  if (!'formula' %in% ...names()) {
    sd_col <- rlang::sym('sd')
  } else {
    sd_col <- rlang::sym(as.character(rlang::list2(...)[['formula']])[2])
  }
  y_tit <- set_y_tit(sd_col)
  gam_reg <- rlang::eval_tidy(
    rlang::call2(fit_gamma_regression, data = data, !!!rlang::dots_list(...))
  )

  c_nms <- sort(unique(data$c))
  clrs  <- viridisLite::turbo(length(c_nms), begin = .9, end = 0)
  clrs2 <- adjustcolor(clrs, alpha.f = .25)

  plt <- data %>%
    tidyr::drop_na(sd) %>%
    ggplot2::ggplot(ggplot2::aes(mean, !!sd_col, color = c)) +
    ggplot2::geom_point(size = 1 / 10)
  plt   <- add_gglayer(plt, c_nms = c_nms, clrs = clrs, gam_reg = gam_reg)
  plt <- plt +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("After Partitioning") +
    ggplot2::scale_color_manual('Trend', values = clrs2) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2))) +
    ggplot2::labs(
      x = expression(bold(bar(y))), y = y_tit
    ) +
    ggplot2::theme(
      legend.position = c(.9, .9)
    )
  if (!is.null(score)) {
    plt +
      ggplot2::annotate(geom = 'text', label = paste0("Score: ", round(score, 5)), x = Inf, y = Inf, hjust = 1, vjust = 1)
  } else {
    plt
  }
}

set_y_tit <- function(sd_col) {
  if (sd_col == 'sd') {
    expression(bold(s))
  } else{
    expression(bold(s[p]))
  }
}

add_gglayer <- function(p, c_nms, clrs, gam_reg) {
  lines <- purrr::map(seq_along(c_nms),
                      function(y) ggplot2::stat_function(
                        fun = ~ stats::predict.glm(
                          gam_reg,
                          newdata = data.frame(mean = .x, c = c_nms[y]),
                          type = "response"
                        ),
                        color = clrs[y])

  )
  p + lines
}
