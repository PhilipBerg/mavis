utils::globalVariables(
  c("alpha", "intu", "intl", "betau", "betal", "res", "everything")
)
#' Mean-Variance Trend Partitioning
#'
#' @description Partitions the data into two mixtures.
#'
#' @param data A `tibble` or `data.frame` to partition
#' @param design_matrix A design matrix for the data (see example)
#' @param formula Formula for the Gamma regression
#' @param eps Size of the integration window
#' @param n Number of integration windows
#' @param verbose If the number of points moved should be output
#' @param eps_scale If the windows should be linear in log-space ("log-linear")
#' or on the real line ("linear").
#'
#' @return A `tibble` or `data.frame` the partitioning vector `c`
#' @export
#'
#' @examples # Please see vignette('baldur') for examples
trend_partitioning <- function(data, design_matrix, formula = sd ~ mean * c, eps = .1, eps_scale = "log-linear", verbose = T, max_iter = 15)
  {
  sd_col <- as.character(formula)[2]
  if (length(eps) == 2) {
    pseq <- purrr::partial(seq, length.out = sum(!is.na(data[[sd_col]])))
    if (eps_scale != 'linear') {
      eps <- log(eps)
      eps <- pseq(eps[1], eps[2])
      eps <- exp(eps)
    } else{
      eps <- pseq(eps[1], eps[2])
    }
    cur_data <- data %>%
      tidyr::drop_na(!!sd_col) %>%
      dplyr::arrange(dplyr::desc(mean))
  } else {
    cur_data <- data
  }
  cur_data <- cur_data %>%
    prep_data_for_clustering(formula, eps = eps) %>%
    run_procedure(formula, eps = eps)
  ll_val <- -Inf
  counter <- 0
  while (cur_data$i > 1) {
    if (verbose) {
      cli::cli_alert_info("Moved {cur_data$i} points")
    }
    counter <- counter + 1
    cur_data <- cur_data$data %>%
      run_procedure(
        formula,
        eps = eps
      )
    if (length(unique(cur_data$data$c)) == 1) {
      cli::cli_alert_danger("Partition converged to a single cluster\nBreaking and returning single partition")
      cli::cli_alert_warning("Try a different window size; ignore if part of windows size estimation")
      break
    }
    ll_val_cur <- cur_data$data %>%
      fit_gamma_regression(formula) %>%
      estimate_gamma_hyperparameters(cur_data$data) %$%
      ll(sd, alpha[1], beta)
    if (ll_val_cur > ll_val) {
      ll_val <- ll_val_cur
    } else if (ll_val_cur < ll_val & counter >= max_iter) {
      if (verbose) print("Global optima increased; breaking")
      break
    }
  }
  if (anyNA(data[[sd_col]])) {
    out <- cur_data$data %>%
      cluster_missing_sd(
        data, design_matrix,
        formula, max(eps)
      )
  }
  else {
    out <- cur_data$data
  }
  out %>%
    dplyr::select(-c(betal, betau, intl, intu, alpha))
}

prep_data_for_clustering <- function(data, formula, eps = .1){
  fc <- as.character(formula)
  data_ms <- data %>%
    tidyr::drop_na(fc[2])
  formula <- stats::as.formula(paste0(fc[2], fc[1], 'mean'))
  gam_reg <- baldur::fit_gamma_regression(data_ms, formula)
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
      c = dplyr::if_else(intl < intu, 'U', 'L')
    )
  i <- sum(tmp != data$c)
  return(
    list(
      data = data,
      i = i
    )
  )
}

resetimate_gamma_pars <- function (data, formula, gm = NULL)
{
  if (is.null(gm)) {
    gm <- baldur::fit_gamma_regression(data, formula)
  }
  data %>% dplyr::mutate(
    alpha = 1/summary(gm)$dispersion,
    betal = estimate_beta(gm, mean, "L", alpha),
    betau = estimate_beta(gm, mean, "U", alpha)
  )
}

add_integrals <-function (data, eps, sd_col)
{
  data %>%
    dplyr::mutate(
      intu = num_int_trapz(!!sd_col, alpha, betau, eps),
      intl = num_int_trapz(!!sd_col, alpha, betal, eps)
    )
}

run_procedure <- function (data, formula, eps)
{
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(
      eps = eps,
      rlang::sym(as.character(formula)[2])
    ) %>%
    clust_itt_norm()
}

num_int_trapz <- function (sd, alpha, beta, eps)
{
  if (length(eps) == 1) {
    eps <- rep(eps, times = length(sd))
  }

  rng <- purrr::map2(sd, eps, ~ c(-.y, .y) + .x)

  purrr::pmap_dbl(list(rng, alpha, beta), ~ pgamma(..1[2], ..2, ..3) - pgamma(..1[1], ..2, ..3))


  # x <- purrr::map(rng, ~ seq(.x[1], .x[2], length.out = n))
  # x <- purrr::map(x, ~ .x[.x>=0])
  # dx <- purrr::map(x, diff)
  #
  # y <- purrr::pmap(list(x, alpha, beta), stats::dgamma)
  # y <- purrr::map(y, ~ (.x[1:(length(.x)-1)] + .x[2:length(.x)]))
  #
  # purrr::map2_dbl(y, dx,
  #                 ~ sum(.x*.y)
  # )
}

cluster_missing_sd <- function(clustered_data, data, design_matrix, formula, eps)
  {
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
    add_integrals(eps, sd_col = sd_col) %>%
    dplyr::mutate(
      c = dplyr::if_else(intl < intu, 'U', 'L'),
      !!sd_col := NA
    ) %>%
    dplyr::bind_rows(clustered_data, .)
}

ll <- function(x, alpha, beta)
  {
  alpha*sum(log(beta)) -
    length(x)*log(gamma(alpha)) +
    (alpha - 1)*sum(log(x)) -
    sum(beta*x)
}
