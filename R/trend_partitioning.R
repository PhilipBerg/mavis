utils::globalVariables(
  c("alpha", "intu", "intl", "betau", "betal", "res", "everything", "s",
    "score", "proposed_score", "tmp", "tmp_global")
)
#' Mean-Variance Trend Partitioning
#'
#' @description Partitions the data into two mixtures.
#'
#' @param data A `tibble` or `data.frame` to partition
#' @param design_matrix A design matrix for the data (see example)
#' @param formula Formula for the Gamma regression
#' @param h Size of the integration window
#' @param verbose If the number of points moved should be output
#'
#' @return A `tibble` or `data.frame` the partitioning vector `c`
#' @export
#'
#' @examples
#' # Normalize data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#'   psrn("identifier")
#' yeast_norm %>%
#'   calculate_mean_sd_trends(design) %>%
#'   trend_partitioning(design)
trend_partitioning <- function(data, design_matrix, formula = sd ~ mean + c, h = .1, verbose = T)
  {
  sd_col <- as.character(formula)[2]
  h <- produce_h_vector(h, data)

  h2_full <- produce_h_vector(h, data)
  na_idx <- is.na(data[[sd_col]])
  h <- h2_full[!na_idx]

  cur_data <- data %>%
    tibble::rowid_to_column("tmp") %>%
    dplyr::filter(!na_idx) %>%
    dplyr::arrange(dplyr::desc(mean))

  cur_data <- cur_data %>%
    prep_data_for_clustering(formula, h = h) %>%
    run_procedure(formula, h = h)

  tracker <- list()
  counter <- 1
  tracker$c[[1]]  <- cur_data$data$c
  tracker$ll[[1]] <- calc_ll(cur_data$data, formula)

  while (cur_data$i > 1) {
    if (verbose) {
      cli::cli_alert_info("Moved {cur_data$i} points")
    }
    if (length(unique(cur_data$data$c)) == 1) {
      cli::cli_alert_danger("Partition converged to a single cluster\nBreaking and returning single partition")
      break
    }
    counter <- counter + 1
    cur_data <- cur_data$data %>%
      run_procedure(
        formula,
        h = h
      )
    tracker$c[[counter]]  <- cur_data$data$c
    tracker$ll[[counter]] <- calc_ll(cur_data$data, formula)
    check_dup <- c(which(duplicated(tracker$c, fromLast = TRUE)), which(duplicated(tracker$c)))
    if (!rlang::is_empty(check_dup) & cur_data$i > 0) {
      mx_ll <- which.max(tracker$ll[check_dup[1]:check_dup[2]])
      cur_data$data$c <- tracker$c[[mx_ll]]
      break
    }
  }
  if (anyNA(data[[sd_col]])) {
    out <- cur_data$data %>%
      cluster_missing_sd(
        data, design_matrix,
        formula, h2_full[na_idx]
      )
  }
  else {
    out <- cur_data$data
  }
  out %>%
    dplyr::arrange(tmp) %>%
    dplyr::select(-c(betal, betau, intl, intu, alpha, tmp))
}

#' Multiple Trend Partitioning
#'
#' Performs multiple trend partitioning on the mean-variance trend and adds a
#' new column classifying each peptide into a cluster.
#'
#'
#' @inheritParams trend_partitioning
#' @importFrom cli cli_alert_danger
#' @importFrom cli cli_alert_warning
#' @importFrom cli cli_alert_success
#' @importFrom cli cli_alert_info
#' @importFrom tidyselect any_of
#' @importFrom dplyr select
#' @importFrom forcats fct_collapse
#'
#' @param verbose Should the procedure be printed (1 for multi partitioning, 2
#' for multi and trend partitioning, 0 for no).
#' @param penalty Penalty function for the penalized likelihood, defaults to [p]
#' @param plot Should a plot of the final solution be printed?
#'
#' @return a `tibble` or `data.frame` with \eqn{\boldsymbol{\hat{C}}} added as a new column.
#' @export
#'
#' @examples
#' # Normalize data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#'   psrn("identifier")
#' yeast_norm %>%
#'   calculate_mean_sd_trends(design) %>%
#'   multi_trend_partitioning(design)
multi_trend_partitioning <- function(data, design_matrix,
                                     formula = sd ~ mean + c, h = 0.1,
                                     verbose = 0 , penalty = p,
                                     plot = TRUE
) {

  verbose_local <- verbose > 1
  verbose <- verbose > 0
  env <- environment()

  data <- data %>%
    tibble::rowid_to_column(var = "tmp_global")
  data_full <- data
  sd_col <- as.character(formula)[2]
  h2_full <- produce_h_vector(h, data)
  na_idx <- is.na(data[[sd_col]])
  h2 <- h2_full[!na_idx]

  data <- data %>%
    dplyr::filter(!na_idx)

  if (!"c" %in% names(data)) {

    best_score <- data %>%
      calculate_objective(sd ~ mean, 1, penalty)

    if (verbose) {
      cli::cli_alert_info(
        "Score single trend:  {best_score}"
      )
      cli::cli_alert_info("Splitting Clusters")
    }

    split_data <- trend_partitioning(data, design_matrix, formula, h, verbose_local)
    update(split_data, env)
  }

  if (best_score > proposed_score) {
    if (verbose) cli::cli_alert_success("No new maxima; breaking")
    if(plot) {
      plt <- plot_gamma(data)
      print(plt)
    }
    return(
      list(
        data  = data,
        score = best_score
      )
    )
  }

  old_best <- -Inf
  while (best_score != old_best) {

    old_best <- best_score

    if (verbose) cli::cli_alert_info("Splitting clusters")
    split_data <- data_splitter(current_data, design_matrix, formula, h, verbose_local)
    update(split_data, env)

    if (verbose) cli::cli_alert_info("Mergin clusters")
    merge_routine()
    if (verbose) cli::cli_alert_info("Global Update")
    gcr <- fit_gamma_regression(current_data, sd ~ mean + c)
    unc_c <- unique(current_data$c)
    split_data_update <- current_data %>%
      estimate_gamma_hyperparameters(gcr, ., design_matrix) %>%
      dplyr::mutate(
        c = purrr::map(unc_c, ~ estimate_beta(gcr, mean, .x, alpha = alpha)) %>%
          list_to_beta(h2, sd, alpha),
        c = unc_c[c]
      )
    update(split_data_update, env)
  }

  if (verbose) cli::cli_alert_success("No new maxima; breaking")

  if (any(na_idx)) {
    gcr <- fit_gamma_regression(current_data, sd ~ mean + c)
    sd_all <- data_full %>%
      select(matches(get_conditions(design_matrix))) %>%
      unlist() %>%
      sd(na.rm = TRUE)
    unc_c <- unique(current_data$c)
    current_data <- data_full %>%
      dplyr::filter(na_idx) %>%
      dplyr::mutate(
        !!sd_col := sd_all
      ) %>%
      dplyr::mutate(
        alpha = (1/summary(gcr)$dispersion),
        c = purrr::map(unc_c, ~ estimate_beta(gcr, mean, .x, alpha = alpha)) %>%
          list_to_beta(h2_full[na_idx], sd, alpha),
        c = unc_c[c],
        !!sd_col := NA
      ) %>%
      dplyr::bind_rows(current_data) %>%
      dplyr::arrange(tmp_global) %>%
      dplyr::select(-tmp_global)
  }

  if(plot) {
    plt <- plot_gamma_partition(current_data, design_matrix, best_score, formula = formula)
    print(plt)
  }
  current_data <- current_data %>%
    select(-any_of(c("alpha", "beta")))
  list(
    data  = current_data,
    score = best_score
  )
}

#' Grid Search for Integration Window
#'
#' Performs a grid search over a set of values of h.
#'
#' @inheritParams multi_trend_partitioning
#'
#' @importFrom dplyr desc
#'
#' @param workers Number of parallel workers to use
#' @param n_h1,n_h2 Number of window sizes to use for `h1`/`h2` if default
#' grid is used.
#' @param h1,h2 A sequence of windows sizes to use in the grid search
#'
#' @return A tibble of integration windows sorted from highest (best) to lowest
#' (worst) penalized likelihood values.
#' @export
#'
#' @examples
#' # Normalize data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#'   psrn("identifier")
#' yeast_norm %>%
#'   calculate_mean_sd_trends(design) %>%
#'   grid_search(design, n_h1 = 2, n_h2 = 2)
grid_search <- function(
    data, design_matrix, penalty = p, workers = 1,
    n_h1 = 5, n_h2 = 5,
    h1 = seq(1e-3, .25, length.out = n_h1)*diff(range(data$sd, na.rm = TRUE)),
    h2 = seq(1e-3, .25, length.out = n_h1)*diff(range(data$sd, na.rm = TRUE)),
    formula = c(sd ~ mean + c),
    plot = TRUE
) {

  out <- tidyr::expand_grid(
    h1, h2, formula
  )

  if(workers > 1) {
    cl <- multidplyr::new_cluster(workers)
    multidplyr::cluster_library(cl,
                                c(
                                  "magrittr",
                                  "mavis",
                                  "purrr",
                                  "dplyr",
                                  "stringr"
                                )
    )
    p <- penalty
    multidplyr::cluster_copy(cl,
                             c(
                               "multi_trend_partitioning",
                               "calculate_objective",
                               "calcuate_objective",
                               "obj",
                               "ll",
                               "p",
                               "data_splitter",
                               "reannotate_c",
                               "data",
                               "design_matrix",
                               "merge_clusters",
                               "make_merge_candidates",
                               "calculate_merge_objective_helper",
                               "reannotate_c_helper",
                               "merge_routine",
                               "update",
                               "tc_tp",
                               "tc_tp_helper",
                               "list_to_beta",
                               "produce_h_vector"
                             )
    )
    out <- out %>%
      multidplyr::partition(cl)
  }
  out <- out %>%
    dplyr::mutate(
      s = purrr::pmap(list(formula, h1, h2),
               ~ multi_trend_partitioning(data, design_matrix, ..1, c(..2, ..3), verbose = 0, penalty = p, plot = F)
      ),
      clustered_data = purrr::map(s, magrittr::use_series, data),
      s = purrr::map_dbl(s, magrittr::use_series, score)
    )

  if (workers > 1) {
    out <- dplyr::collect(out)
    rm(cl)
    gc()
  }
  out <- out %>%
    dplyr::arrange(desc(s))
  if(plot) {
    plt <- out %>%
      ggplot2::ggplot(ggplot2::aes(h1, h2, fill = s)) +
      ggplot2::geom_raster(interpolate = T) +
      ggplot2::facet_grid(. ~ as.character(formula)) +
      ggplot2::scale_fill_viridis_c(option = "H") +
      ggplot2::coord_cartesian(
        xlim = range(h1), ylim = range(h2)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
    print(plt)
  }
  out
}

#' Penalty function
#'
#' Used to calculate and subtract a penalty from the likelihood.
#'
#' @param k The number of clusters
#' @param n The number of data points
#'
#' @return A real number of the penalty of the likelihood
#' @export
#'
#' @examples
#' p(1, 10)
p <- function(k, n) {
  sqrt(n)*exp(k)
}

produce_h_vector <- function(h, data) {
  leng_h <- c(length(h), nrow(data))
  if (leng_h[1] != leng_h[2]) {
    if (leng_h[1] == 1) {
      h <- rep(h, leng_h[2])
    } else{
      h <- seq(h[1], h[2], length.out = leng_h[2])
    }
  }
  h
}

prep_data_for_clustering <- function(data, formula, h = .1){
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

resetimate_gamma_pars <- function(data, formula, gm = NULL) {
  if (is.null(gm)) {
    gm <- baldur::fit_gamma_regression(data, formula)
  }
  data %>% dplyr::mutate(
    alpha = 1/summary(gm)$dispersion,
    betal = estimate_beta(gm, mean, "L", alpha),
    betau = estimate_beta(gm, mean, "U", alpha)
  )
}

add_integrals <- function(data, h, sd_col) {
  data %>%
    dplyr::mutate(
      intu = num_int_trapz(!!sd_col, alpha, betau, h),
      intl = num_int_trapz(!!sd_col, alpha, betal, h)
    )
}

run_procedure <- function(data, formula, h) {
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(
      h = h,
      rlang::sym(as.character(formula)[2])
    ) %>%
    clust_itt_norm()
}

num_int_trapz <- function(sd, alpha, beta, h) {
  if (length(h) == 1) {
    h <- rep(h, times = length(sd))
  }

  rng <- purrr::map2(sd, h, ~ c(-.y, .y) + .x)

  purrr::pmap_dbl(list(rng, alpha, beta), ~ pgamma(..1[2], ..2, ..3) - pgamma(..1[1], ..2, ..3))

}

cluster_missing_sd <- function(
    clustered_data,
    data,
    design_matrix,
    formula,
    h
) {
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
    add_integrals(h, sd_col = sd_col) %>%
    dplyr::mutate(
      c = dplyr::if_else(intl < intu, 'U', 'L'),
      !!sd_col := NA
    ) %>%
    dplyr::bind_rows(clustered_data, .)
}

ll <- function(x, alpha, beta) {
  alpha*sum(log(beta)) -
    length(x)*log(gamma(alpha)) +
    (alpha - 1)*sum(log(x)) -
    sum(beta*x)
}

calc_ll <- function(data, formula) {
  data %>%
    fit_gamma_regression(formula) %>%
    estimate_gamma_hyperparameters(data) %$%
    ll(sd, alpha[1], beta)
}

merge_routine <- function() {
  x <- parent.frame()
  while (T) {
    merge_list   <- make_merge_candidates(x$split_data$c)
    merge_scores <- purrr::map_dbl(merge_list, calculate_merge_objective_helper, x$split_data, x$formula, x$penalty)
    idx          <- which.max(merge_scores)
    x$split_data <- merge_clusters(x$split_data, merge_list[[idx]]) %>%
      reannotate_c()
    update(x$split_data, x)
    if (length(merge_list) == 1) break
  }
}

make_merge_candidates <- function(clusters) {
  sc <- clusters %>%
    unique() %>%
    sort()
  purrr::map(seq_len(length.out = length(sc) - 1), ~ c(sc[.x], sc[.x + 1]))
}

data_splitter <- function(data, design_matrix, formula = sd ~ mean + c,
                          h = 0.1, verbose = T, max_iter) {
  data %>%
    split.data.frame(.$c) %>%
    purrr::map(
      tc_tp, design_matrix, formula, h, verbose, max_iter = max_iter
    ) %>%
    purrr::imap(
      ~ dplyr::mutate(.x, c = stringr::str_c(.y, c, sep = "."))
    ) %>%
    dplyr::bind_rows() %>%
    reannotate_c()
}

reannotate_c <- function(data) {
  data$c <- reannotate_c_helper(data$c)
  data
}

reannotate_c_helper <- function(c) {
  replacer <- c %>%
    unique() %>%
    sort() %>%
    setNames(letters[1:length(.)], .)
  stringr::str_replace_all(c, replacer)
}

list_to_beta <- function(lst, h, sd, alpha) {
  purrr::map(lst, ~ num_int_trapz(sd, alpha, .x, h)) %>%
    do.call(cbind, .) %>%
    apply(1, which.max)
}

calcuate_objective <- function(data, formula, penalty){
  k <- data$c %>%
    unique() %>%
    length()
  if (k > 1) {
    gr <- data %>%
      fit_gamma_regression(formula)

  } else {
    gr <- data %>%
      fit_gamma_regression(sd ~ mean)
  }
  gr %>%
    estimate_gamma_hyperparameters(data) %$%
    obj(sd, alpha, beta, k, penalty)
}

calculate_objective <- function(data, formula, k, penalty) {
  data %>%
    fit_gamma_regression(formula) %>%
    estimate_gamma_hyperparameters(data) %$%
    obj(sd, alpha, beta, k, penalty)
}


tc_tp <- function(data, design_matrix, formula = sd ~ mean + c, h = 0.1,
                  verbose = T, max_iter) {
  tryCatch(
    trend_partitioning(
      data, design_matrix, formula = formula, h = h, verbose = verbose
    ),
    error = \(x) tc_tp_helper(data), finally = function(e) e
  )
}

tc_tp_helper <- function(data) {
  cli::cli_alert_danger("Partition converged to a single cluster\nBreaking and returning single partition")
  cli::cli_alert_warning("Try a different window size; ignore if part of windows size estimation")
  data
}

obj <- function(x, alpha, beta, k, penalty = p) {
  alpha <- alpha[1]
  n <- length(x)
  ll(x, alpha, beta) - penalty(k, n)
}

merge_clusters <- function(data, who_to_merge) {
  data %>%
    dplyr::mutate(
      c = forcats::fct_collapse(c, !!who_to_merge[1] := who_to_merge)
    ) %>%
    reannotate_c()
}

make_merge_candidates <- function(clusters) {
  sc <- clusters %>%
    unique() %>%
    sort()
  purrr::map(seq_len(length.out = length(sc) - 1), ~ c(sc[.x], sc[.x + 1]))
}

calculate_merge_objective_helper <- function(who_to_merge, data, formula, penalty) {
  data <- data %>%
    merge_clusters(who_to_merge = who_to_merge)
  calcuate_objective(data, formula, penalty)
}

update <- function(data, x) {

  x$proposed_score <- data %>%
    calcuate_objective(x$formula, x$penalty)

  if (x$proposed_score >= x$best_score) {

    if (x$verbose) {
      cli::cli_alert_info(
        "New best score: {x$proposed_score}"
      )
      cli::cli_alert("Current number of clusters: {length(unique(data$c))}")
    }
    x$best_score <- x$proposed_score
    x$current_data  <- data

  }
}

