#' Imputation
#' @param data a `data.frame` to perform the imputation on, missing values should
#' be `NA`.
#' @param design a design or model matrix as produced by
#'  \code{\link[stats]{model.matrix}} with column names corresponding to the
#' @param formula Regression formula for the pooled variance and mean trend
#' @param workers Number of parallel workers to use when building the forest
#' @param imp_pars Imputation parameters to use; as inferred by [get_imputation_pars].
#' @param ... Additional arguments to [ranger::ranger].
#' @name imp
NULL

utils::globalVariables(c("matches", "sd_p", "predictions"))

#' Estimates the parameters for the imputation model.
#' @rdname imp
#'
#' @return a `list` containing the imputation parameters and corresponding index.
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important.
#' colnames(design) <- paste0("ng", c(50, 100))
#' \donttest{
#' # Estimate the imputation parameters
#' imp_pars <- yeast %>%
#'   psrn("identifier") %>%
#'   get_imputation_pars(design)
#' }
get_imputation_pars <- function(data, design, formula = sd_p ~ mean, workers = 1, ...) {
  cat("Estimating Imputation Paramters\n")
  condi <- get_conditions(design)
  mis_vals <- data %>%
    dplyr::select(matches(condi)) %>%
    purrr::map(~which(is.na(.x)))

  mis_vals <- mis_vals[purrr::map_lgl(mis_vals, ~ length(.x) != 0)]

  miss_count <- data %>%
    dplyr::select(matches(condi)) %>%
    split.default(
      stringr::str_extract(colnames(.), condi)
    ) %>%
    purrr::map(
      ~ is.na(.x)
    ) %>%
    purrr::map(
      rowSums
    ) %>%
    setNames(
      paste0('misscount_', seq_along(.))
    )

  sd_pooled <- data %>%
    calculate_mean_sd_trends(design, "pooled")
  gam_reg <- fit_gamma_regression(sd_pooled, formula)
  sd_pooled <- sd_pooled %>%
    magrittr::use_series(sd_p) %>%
    sqrt()

  imp_order <- purrr::map_dbl(mis_vals, length) %>%
    sort()
  mis_vals <- mis_vals[names(imp_order)]

  if (workers != 1) {
    cl <- parallel::makeCluster(
      min(sum(design) - 1, workers)
    )
    doParallel::registerDoParallel(cl)
  }
  dif <- numeric()
  n_diff <- Inf

  imp_mat <- data %>%
    dplyr::select(matches(condi))


  mean_vals <- imp_mat %>%
    purrr::map_dbl(mean, na.rm = T)
  for (i in names(mis_vals)) {
    data[mis_vals[[i]],i] <- imp_mat[mis_vals[[i]],i] <- mean_vals[i]
  }

  c_check <- stringr::str_detect(
    as.character(formula), 'c'
  ) %>%
    any()
  if (c_check) {
    c_vals <- purrr::map(mis_vals, ~ data$c[.x])
  }
  rf_pars <- list(mtry = floor(sum(design)/3))
  in_pars <- rlang::list2(...)
  rf_pars[names(in_pars)] <- in_pars[names(in_pars)]
  rf <- purrr::partial(ranger::ranger, !!!rf_pars)

  while (TRUE) {
    tic <- Sys.time()
    den <- num <- 0
    for(i in names(mis_vals)) {
      form <- formula(
        paste0(i, '~ .')
      )
      idx <- mis_vals[[i]]

      pred <- rf(
        form, dplyr::slice(imp_mat, -idx),
        num.threads = workers
      ) %>%
        stats::predict(dplyr::slice(imp_mat, idx)) %>%
        magrittr::use_series(predictions)

      den <- den + sum(
        (imp_mat[idx,i] - pred)^2
      )

      num <- num + sum(pred^2)

      data[idx,i] <- imp_mat[idx,i] <- pred
    }
    dif <- den / num
    if (n_diff > dif) {
      cat(
        'Previous error:', n_diff, '\t>\tCurrent error:', sum(dif), '\n'
      )

      n_diff <- dif
      imp_out <- data
      print_it_time(tic)
    } else {
      cat(
        'Previous error:', n_diff, '\t<\tCurrent error:', sum(dif), '\tBreaking \n'
      )
      sd_error <- means <- list()
      for (i in names(mis_vals)) {
        idx <- mis_vals[[i]]
        means[[i]] <- imp_mat[[i]][idx]
        new_data <- data.frame(mean = means[[i]])
        if (c_check) {
          new_data <- dplyr::bind_cols(new_data, c = c_vals[[i]])
        }
        trend_sd <- stats::predict.glm(
          gam_reg,
          new_data,
          type = 'response'
        ) %>%
          sqrt()
        sd_error[[i]] <- dplyr::if_else(
          !is.na(sd_pooled[idx]),
          sd_pooled[idx] * trend_sd,
          trend_sd^2
        )
      }
      print_it_time(tic)
      break
    }
  }
  if (workers != 1) {
    parallel::stopCluster(cl)
    rm(cl)
    gc()
  }
  return(
    list(
      mis_vals = mis_vals,
      means = means,
      sd_error = sd_error
    )
  )
}

#' Performs a single imputation run and returns the data with NA values replaced
#' by imputed values.
#'
#' @rdname imp
#'
#' @export
#' @importFrom ranger ranger
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @return a `data.frame` with `NA` values replaced by imputed values.
#'
#' @examples
#' \donttest{
#' # Run the imputation
#' yeast %>%
#'   single_imputation(design, imp_pars = imp_pars)
#' }
single_imputation <- function(data,
                              design,
                              formula = sd_p ~ mean,
                              workers = 1,
                              imp_pars = NULL,
                              ...) {
  if (is.null(imp_pars)) {
    imp_pars <- get_imputation_pars(data, design, formula, workers, ...)
  }
  for (i in names(imp_pars$mis_vals)) {
    data[imp_pars$mis_vals[[i]],i] <- imp_pars$means[[i]]
  }
  return(data)
}

impute <- function(data, imp_pars) {
  for (i in names(imp_pars$mis_vals)) {
    data[imp_pars$mis_vals[[i]],i] <- stats::rnorm(
      n    = length(imp_pars$mis_vals[[i]]),
      mean = imp_pars$means[[i]],
      sd   = imp_pars$sd_error[[i]]
    )
  }
  return(data)
}

print_it_time <- function(tic) {
  toc <- Sys.time() - tic
  cat(
    "Iteration time:\n",
    format(toc), '\n\n'
  )
}
