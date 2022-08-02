utils::globalVariables(
  c("alpha", "intu", "intl", "betau", "betal", "res", "everything")
)
#' Mean-Variance Trend Partitioning
#'
#' @description Partitions the data into two mixtures.
#' @param data A `tibble` or `data.frame` to partition
#' @param design_matrix A design matrix for the data (see example)
#' @param formula Formula for the Gamma regression
#' @param eps Size of the integration window
#' @param n Number of integration windows
#' @param verbose If the number of points moved should be output
#'
#' @return A `tibble` or `data.frame` the partitioning vector `c`
#' @export
#'
#' @examples # Please see vignette('baldur') for examples
trend_partitioning <- function(data, design_matrix, formula = sd ~ mean + c, eps = .1, n = 1000, verbose = T){
  cur_data <- data %>%
    prep_data_for_clustering(design_matrix, eps = eps, n = n) %>%
    run_procedure(formula, eps = eps, n = n)
  while (cur_data$i > 1) {
    if(verbose){
      print(
        paste('Moved', cur_data$i, 'points')
      )
    }
    cur_data <- cur_data$data %>%
      run_procedure(formula, eps = eps, n = n)
  }
  if(any(is.na(cur_data$data$sd))) {
  cur_data$data %>%
    cluster_missing_sd(data, design_matrix, formula, eps, n) %>%
    dplyr::select(-c(betal, betau, intl, intu, alpha))
  } else {
    cur_data$data
  }
}

prep_data_for_clustering <- function(data, design_matrix, eps = .1, n = 1000){
  data_ms <- data %>%
    tidyr::drop_na(sd)
  gam_reg <-  stats::glm(sd ~ mean, stats::Gamma(log), data_ms)
  data_ms %>%
    dplyr::mutate(
      res = stats::residuals(gam_reg),
      c = dplyr::if_else(res<0, 'L', 'U')
    ) %>%
    dplyr::select(-res)
}

clust_itt_norm <- function(data, reg){
  tmp <- data$c
  data <- data %>%
    dplyr::mutate(
      c = dplyr::if_else(intl<intu, 'U', 'L')
    )
  i <- sum(tmp != data$c)
  return(
    list(
      data = data,
      i = i
    )
  )
}

resetimate_gamma_pars <- function(data, formula, gm = NULL){
  if(is.null(gm)){
    gm <- stats::glm(formula, stats::Gamma(log), data)
  }
  data %>%
    dplyr::mutate(
      alpha = (1/summary(gm)$dispersion),
      betal = estimate_beta(gm, mean, 'L', alpha),
      betau = estimate_beta(gm, mean, 'U', alpha)
    )
}

add_integrals <- function(data, eps = .1, n = 1000){
  data %>%
    dplyr::mutate(
      intu = purrr::pmap_dbl(list(sd, alpha, betau), num_int_trapz, eps, n),
      intl = purrr::pmap_dbl(list(sd, alpha, betal), num_int_trapz, eps, n)
    )
}

run_procedure <- function(data, formula, eps = .1, n = 1000){
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(eps = eps, n = n) %>%
    clust_itt_norm()
}

num_int_trapz <- function(sd, alpha, beta, eps, n){
  rng <- c(-eps, eps) + sd
  x <- seq(rng[1], rng[2], length.out = n)
  dx <- diff(x)
  y <- stats::dgamma(x, alpha, beta)
  sum(c(y[1]/2, y[2:(length(y)-2)], y[length(y)-1]/2)*dx)
}

estimate_beta <- function(model, mean, c, alpha){
  alpha / stats::predict.glm(model, newdata = data.frame(mean = mean, c = c), type = 'response')
}

cluster_missing_sd <- function(clustered_data, data, design_matrix, formula, eps, n){
  sigma_all <- data %>%
    dplyr::select(dplyr::matches(get_conditions(design_matrix))) %>%
    unlist() %>%
    sd(na.rm = T)
  gam_mod <- stats::glm(sd ~ mean + c, stats::Gamma(log), clustered_data)
  data %>%
    dplyr::filter(is.na(sd)) %>%
    dplyr::mutate(
      sd = sigma_all
    ) %>%
    resetimate_gamma_pars(formula, gam_mod) %>%
    add_integrals(eps, n) %>%
    dplyr::mutate(
      c = dplyr::if_else(intl<intu, 'U', 'L'),
      sd = NA
    ) %>%
    dplyr::bind_rows(clustered_data)
}
