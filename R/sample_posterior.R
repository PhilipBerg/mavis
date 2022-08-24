utils::globalVariables(c("tmp"))
#' Sample Posterior
#'
#' @description Function to sample the posterior of the Bayesian decision model.
#' Note that, `sample_posterior` will remove missing values and adjust the inputs
#' to the Stan model accordingly. If all values in a condition is missing for a
#' row/peptide, then that will be removed for that row/peptide.
#' Still, if there is only one condition that has observed values, `sample_posterior`
#' will crash. Please make sure to impute or filter out these cases.
#'
#' @param data A `tibble` or `data.frame` with alpha,beta priors annotated
#' @param id_col_name A character of the id column
#' @param design_matrix A design matrix for the data (see example)
#' @param contrast_matrix A contrast matrix of the decisions to test such that the first column is the index of the column in the design matrix that should be compared to the second column
#' @param uncertainty_matrix A matrix of observation uncertainty
#' @param bayesian_model Which Bayesian model to use. Currently only on internal model allowed, this argument is for forward compatibility
#' @param clusters The number of parallel threads to run on.
#' @param warmup Number of warm-up samples to draw from the model
#' @param iter Total number of iterations to draw (i.e., iter - warmup = posterior draws)
#' @param chains The number of chains to sample from
#'
#' @return A `tibble` or `data.frame` annotated  with statistics of the posterior
#' @export
#'
#' @examples # Please see vignette('baldur') for examples
sample_posterior <- function(
    data,
    id_col_name,
    design_matrix,
    contrast_matrix,
    uncertainty_matrix,
    bayesian_model = stanmodels$uncertainty_model,
    clusters = 1,
    warmup = 1000,
    iter = 2000,
    chains = 4) {
  N <- sum(design_matrix)
  K <- ncol(design_matrix)
  C <- nrow(contrast_matrix)
  pmap_columns <- rlang::expr(list(!!rlang::sym(id_col_name), alpha, beta))
  ori_data <- data
  if(clusters != 1) {
    #### ####
    cl <- multidplyr::new_cluster(clusters)
    multidplyr::cluster_library(cl,
                                c(
                                  "dplyr",
                                  "tidyr",
                                  "purrr",
                                  "tibble",
                                  "stringr",
                                  "magrittr",
                                  "rstan",
                                  'StanHeaders'
                                )
    )
    multidplyr::cluster_copy(cl,
                             c(
                               "bayesian_testing",
                               "id_col_name",
                               "design_matrix",
                               "uncertainty_matrix",
                               "generate_stan_data_input",
                               "stan_summary",
                               "pmap_columns",
                               "ori_data",
                               "N",
                               "K",
                               "C",
                               "contrast_matrix",
                               'bayesian_model'
                             )
    )
    model_check <- stringr::word(deparse(substitute(bayesian_model)), 2, sep = '\\$')
    if(!is.na(model_check)) {
      multidplyr::cluster_copy(cl,
                               paste0(
                                 "rstantools_model_",
                                 stringr::word(deparse(substitute(bayesian_model)), 2, sep = '\\$')
                               )
      )
    }
    data <- ori_data %>%
      multidplyr::partition(cl)
    on.exit({rm(cl); gc()})
  }
  data %>%
    dplyr::mutate(
      tmp = purrr::pmap(!!pmap_columns,
                        bayesian_testing,
                        ori_data, id_col_name, design_matrix, contrast_matrix, bayesian_model, N, K, C, uncertainty_matrix,
                        warmup, iter, chains
      )
    ) %>%
    dplyr::collect() %>%
    tidyr::unnest(tmp)
}

generate_stan_data_input <- function(id, id_col_name, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex) {
  row <- data[
    data[[id_col_name]] == id, stringr::str_detect(names(data), condi_regex)
  ]
  xbar <- purrr::map_dbl(purrr::set_names(condi, condi),
                         ~as.numeric(row[stringr::str_which(colnames(row), .x)]) %>%
                           mean(na.rm = T)
  )[condi]
  if(is.null(uncertainty)) {
    u = rep(1,N)
  }else{
    u = uncertainty[id,]
  }
  list(
    N = N,
    K = K,
    C = C,
    x = design_matrix,
    y = as.numeric(row),
    c = comparison,
    alpha = alpha,
    beta_gamma = beta,
    xbar = xbar,
    u = u
  )
}

modify_stan_input_for_missingingess <- function(stan_input) {
  na_idx <- which(is.na(stan_input$y))
  stan_input$x <- stan_input$x[-na_idx,]
  stan_input$y <- stan_input$y[-na_idx]
  stan_input$u <- stan_input$u[-na_idx]
  stan_input$N <- stan_input$N - length(na_idx)
  miss_condi_idx <- which(colSums(stan_input$x) == 0)
  lng_miss <- length(miss_condi_idx)
  if(lng_miss != 0) {
    stan_input$x <- stan_input$x[, -miss_condi_idx, drop = F]
    stan_input$xbar <- stan_input$xbar[-miss_condi_idx]
    stan_input$K <- stan_input$K - lng_miss
    stan_input$c <- stan_input$c[
      -which(stan_input$c == miss_condi_idx, arr.ind = T)[,1],, drop = F
    ]
    for (i in miss_condi_idx) {
      stan_input$c[stan_input$c > i] <- stan_input$c[stan_input$c > i] - 1
    }
    stan_input$C <- nrow(stan_input$c)
  }
  return(stan_input)
}

stan_summary <- function(fit, condi, C, dat) {
  lfc_pars <- paste0('y_diff[', seq_len(C), ']')
  err_pars <- paste0('error[', seq_len(C), ']')
  summ <- rstan::summary(fit) %>%
    magrittr::use_series(summary)
  summ <- summ[
    rownames(summ) %in% c(lfc_pars, err_pars, 'sigma', 'lp__'),
    colnames(summ) %in% c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat', '50%')
  ]
  summ <- summ[1:nrow(summ), c(1:2, 4:6, 3)]
  dplyr::tibble(
    comparison = purrr::map2_chr(condi[dat$c[,1]], condi[dat$c[,2]], ~stringr::str_flatten(c(.x, .y), collapse = ' vs ')),
    err = summ[rownames(summ) %in% err_pars, 1],
    err_025 = summ[rownames(summ) %in% err_pars, 2],
    err_975 = summ[rownames(summ) %in% err_pars, 3],
    err_eff = summ[rownames(summ) %in% err_pars, 4],
    err_rhat = summ[rownames(summ) %in% err_pars, 5],
    lfc = summ[rownames(summ) %in% lfc_pars, 1],
    lfc_025 = summ[rownames(summ) %in% lfc_pars, 2],
    lfc_975 = summ[rownames(summ) %in% lfc_pars, 3],
    lfc_eff = summ[rownames(summ) %in% lfc_pars, 4],
    lfc_rhat = summ[rownames(summ) %in% lfc_pars, 5],
    sigma = summ[rownames(summ) == 'sigma', 1],
    sigma_025 = summ[rownames(summ) == 'sigma', 2],
    sigma_975 = summ[rownames(summ) == 'sigma', 3],
    sigma_eff = summ[rownames(summ) == 'sigma', 4],
    sigma_rhat = summ[rownames(summ) == 'sigma', 5],
    lp__ = summ[rownames(summ) == 'lp__', 1],
    lp___025 = summ[rownames(summ) == 'lp__', 2],
    lp___975 = summ[rownames(summ) == 'lp__', 3],
    lp___eff = summ[rownames(summ) == 'lp__', 4],
    lp___rhat = summ[rownames(summ) == 'lp__', 5]
  )
}

bayesian_testing <- function(
    id,
    alpha,
    beta,
    data,
    id_col_name,
    design_matrix,
    comparison,
    model,
    N,
    K,
    C,
    uncertainty = NULL,
    warmup,
    iter,
    chains) {
  condi <- colnames(design_matrix)
  condi_regex <- paste0(condi, collapse = '|')
  dat <- generate_stan_data_input(id, id_col_name, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex)
  if(any(is.na(dat$y))) {
    dat <- modify_stan_input_for_missingingess(dat)
    condi <- colnames(dat$x)
    condi_regex <- paste0(condi, collapse = '|')
  }
  fit <- purrr::quietly(rstan::sampling)(object = model, data = dat, verbose = F, refresh = 0, warmup = warmup, iter = iter, chains = chains)
  if(length(fit$warnings) != 0) {
    cat(fit$warnings, '\n')
    cat('Divergent transition for:\n', id, '\n')
  }
  fit$result %>%
    stan_summary(condi, C, dat)
}
