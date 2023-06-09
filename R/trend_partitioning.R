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
    prep_data_for_clustering(formula, eps = eps, n = n) %>%
    run_procedure(formula, eps = eps, n = n)
  sd_col <- as.character(formula)[2]
  while (cur_data$i > 1) {
    if(verbose){
      print(
        paste('Moved', cur_data$i, 'points')
      )
    }
    cur_data <- cur_data$data %>%
      run_procedure(formula, eps = eps, n = n)
  }
  if(anyNA(data[[sd_col]])) {
    out <- cur_data$data %>%
      cluster_missing_sd(data, design_matrix, formula, eps, n)
  } else {
    out <-  cur_data$data %>%
      dplyr::select(-c(betal, betau, intl, intu, alpha))
  }
  class(out) <- c('gmr', class(out))
}

prep_data_for_clustering <- function(data, formula, eps = .1, n = 1000){
  fc <- as.character(formula)
  data_ms <- data %>%
    tidyr::drop_na(fc[2])
  formula <- as.formula(paste0(fc[2], fc[1], 'mean'))
  gam_reg <-  stats::glm(formula, stats::Gamma(log), data_ms)
  data_ms %>%
    dplyr::mutate(
      res = stats::residuals(gam_reg),
      c = dplyr::if_else(res < 0, 'L', 'U')
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

add_integrals <-function (data, eps, n = 1000, sd_col)
{
  data %>%
    dplyr::mutate(
      intu = num_int_trapz(!!sd_col, alpha, betau, eps, n),
      intl = num_int_trapz(!!sd_col, alpha, betal, eps, n)
    )
}

run_procedure <- function(data, formula, eps = .1, n = 1000){
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(eps = eps, n = n, rlang::sym(as.character(formula)[2])) %>%
    clust_itt_norm()
}

num_int_trapz <- function (sd, alpha, beta, eps, n)
{
  if (length(eps) == 1) {
    eps <- rep(eps, times = length(sd))
  }
  rng <- purrr::map2(sd, eps, ~ c(-.y, .y) + .x)
  x <- purrr::map(rng, ~ seq(.x[1], .x[2], length.out = n))
  dx <- purrr::map(x, diff)
  y <- purrr::pmap(list(x, alpha, beta), stats::dgamma)
  y <- purrr::map(y, ~ (.x[1:(length(.x)-1)] + .x[2:length(.x)]))
  purrr::map2_dbl(y, dx,
                  ~ sum(.x*.y)
  )
}

cluster_missing_sd <- function(clustered_data, data, design_matrix, formula, eps, n){
  sigma_all <- data %>%
    dplyr::select(dplyr::matches(get_conditions(design_matrix))) %>%
    unlist() %>%
    sd(na.rm = T)
  gam_mod <- stats::glm(formula, stats::Gamma(log), clustered_data)
  sd_col <- rlang::sym(as.character(formula)[2])
  data %>%
    dplyr::filter(is.na(!!sd_col)) %>%
    dplyr::mutate(
      !!sd_col := sigma_all
    ) %>%
    resetimate_gamma_pars(formula, gam_mod) %>%
    add_integrals(eps, n, sd_col = sd_col) %>%
    dplyr::mutate(
      c = dplyr::if_else(intl < intu, 'U', 'L'),
      !!sd_col := NA
    ) %>%
    dplyr::bind_rows(clustered_data, .)
}
