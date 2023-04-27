#'@importFrom Rdpack reprompt
utils::globalVariables(c(".", "sd", "model"))
#' Function for plotting the mean-variance gamma regressions
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the precision weights and imputation.
#'
#'
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param ... Additional arguments to \code{\link[mavis]{trend_partitioning}}.
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
#' yeast <- psrn(yeast_prog, "identifier")
#'
#' # Generate the plots
#' \dontrun{
#' plot_gamma_regression(yeast, design, verbose = FALSE)
#' }
plot_gamma_regression <- function(data, design, ...) {
  data <- data %>%
    calculate_mean_sd_trends(design)
  base_plot <- data %>%
    plot_gamma()
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

plot_mean_sd_trend <- function(data) {
  data %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd)) +
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
      x = expression(bold(bar(y))), y = expression(bold(s))
    )
}

#' Function for plotting the gamma regression used for mean-variance normalization
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the precision weights
#'
#' @param data The data to use for producing the plots.
#'
#' @return a plot with the mean-variance trend used for the precision
#' weights.
#' @export
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- psrn(yeast_prog, "identifier")
#'
#' # Generate the plot
#' yeast %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma()
plot_gamma <- function(data) {
  data %>%
    tidyr::drop_na(sd) %>%
    plot_mean_sd_trend() +
    ggplot2::ggtitle("Before Partitioning")
}

#' Function for plotting the gamma regression used for imputation
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends for the imputation.
#'
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param ... Additional arguments to \code{\link[mavis]{trend_partitioning}}.
#'
#' @return a plot with the mean-variance trend used for the trend lines used in
#'  the imputation.
#' @export
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- psrn(yeast_prog, "identifier")
#'
#' # Generate the plot
#' \dontrun{
#' yeast %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma_partition(design, verbose = FALSE)
#' }
plot_gamma_partition <- function(data, design, ...) {
  if(!'c' %in% names(data)){
    data <- data %>%
      trend_partitioning(design, ...)
  }
  trend_colors <- purrr::set_names(viridisLite::turbo(2, end = .75), c('Lower', 'Upper'))
  gam_reg <- rlang::eval_tidy(
    rlang::call2(fit_gamma_regression, data = data, !!!rlang::dots_list(...))
  )
  data %>%
    dplyr::mutate(
      c = plyr::mapvalues(c, c('L', 'U'), c('Lower', 'Upper'))
    ) %>%
    tidyr::drop_na(sd) %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd, color = c)) +
    ggplot2::geom_point(size = 1 / 10) +
    ggplot2::stat_function(
      fun = ~stats::predict.glm(gam_reg, newdata = data.frame(mean = .x, c = 'L'), type = 'response'), color = 'blue'
    ) +
    ggplot2::stat_function(
      fun = ~stats::predict.glm(gam_reg, newdata = data.frame(mean = .x, c = 'U'), type = 'response'), color = 'red'
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("After Partitioning") +
    ggplot2::scale_color_manual('Trend', values = trend_colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2))) +
    ggplot2::labs(
      x = expression(bold(bar(y))), y = expression(bold(s))
    ) +
    ggplot2::theme(
      legend.position = c(.9, .9)
    )
}
