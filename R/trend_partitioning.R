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
trend_partitioning <- function(data, design_matrix, formula = sd ~ mean * c, eps = .1, eps_scale = "log-linear", n = 1000, verbose = T){
  sd_col <- as.character(formula)[2]
  if (length(eps) == 2) {
    lo <- sum(!is.na(data[[sd_col]]))
    if (eps_scale != 'linear') {
      eps <- log(eps)
      eps <- seq(eps[1], eps[2], lo)
      eps <- exp(eps)
    } else{
      eps <- seq(eps[1], eps[2], lo)
    }
    cur_data <- data %>%
      tidyr::drop_na(!!sd_col) %>%
      dplyr::arrange(dplyr::desc(mean))
  } else {
    cur_data <- data
  }
  cur_data <- cur_data %>%
    prep_data_for_clustering(formula, eps = eps, n = n) %>%
    run_procedure(formula, eps = eps, n = n)
  while (cur_data$i > 1) {
    if (verbose) {
      print(paste("Moved", cur_data$i, "points"))
    }
    cur_data <- cur_data$data %>%
      run_procedure(
        formula,
        eps = eps, n = n
      )
  }
  if (anyNA(data[[sd_col]])) {
    out <- cur_data$data %>%
      cluster_missing_sd(
        data, design_matrix,
        formula, max(eps), n
      )
  }
  else {
    out <- cur_data$data
  }
  out <- out %>%
    dplyr::select(-c(betal, betau, intl, intu, alpha))
  class(out) <- c('gmr', class(out))
  return(out)
}

prep_data_for_clustering <- function(data, formula, eps = .1, n = 1000){
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

add_integrals <-function (data, eps, n = 1000, sd_col)
{
  data %>%
    dplyr::mutate(
      intu = num_int_trapz(!!sd_col, alpha, betau, eps, n),
      intl = num_int_trapz(!!sd_col, alpha, betal, eps, n)
    )
}

run_procedure <- function (data, formula, eps, n = 1000)
{
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(
      eps = eps,
      n = n, rlang::sym(as.character(formula)[2])
    ) %>%
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

estimate_beta <- function(model, mean, c, alpha, ...){
  reg_vars <- model %>%
    stats::formula() %>%
    all.vars()
  reg_vars <- reg_vars[!reg_vars %in% c("mean", "sd", "sd_p")]
  if(length(reg_vars) != 0) {
    reg_vars <- reg_vars %>%
      paste0(., " = ", ., collapse = ", ") %>%
      paste0(', ', .)
  }
  reg_vars <- paste0("newdata = data.frame(mean = mean", reg_vars, ")")
  nd <- rlang::parse_expr(reg_vars)
  alpha / rlang::eval_tidy(
    rlang::call2(
      stats::predict.glm,
      object = model,
      newdata = nd,
      type = 'response'
    )
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
