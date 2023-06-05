utils::globalVariables("mean_condi")
#' Single imputation
#'
#' Performs a single imputation run and returns the data with NA values replaced
#' by imputed values.
#'
#' @param data a `data.frame` to perform the imputation on, missing values should
#' be `NA`.
#' @param design a design or model matrix as produced by
#'  \code{\link[stats]{model.matrix}} with column names corresponding to the
#'  different conditions.
#' @param gam_reg a gamma regression model
#' @return a `data.frame` with `NA` values replaced by imputed values.
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_weights
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast <- yeast_prog %>%
#'   # Normalize and log-transform the data
#'   psrn("identifier")
#'
#' # Run the imputation
#' \dontrun{
#' yeast %>%
#'   single_imputation(design)
#' # Run with partitioning
#' }
single_imputation <- function(data,
                              design,
                              formula = sd_p ~ mean,
                              workers = 1) {
  imp_pars <- get_imputation_pars(data, design, formula, workers)
  impute(data, imp_pars)
}

impute <- function(data, imp_pars) {
  for (i in names(imp_pars$mis_vals)) {
    data[imp_pars$mis_vals[[i]],i] <- rnorm(
      n    = length(imp_pars$mis_vals[[i]]),
      mean = imp_pars$means[[i]],
      sd   = imp_pars$sd_error[[i]]
    )
  }
  return(data)
}


get_imputation_pars <- function(data, design, formula = sd_p ~ mean, workers = 1) {
  cat("Estimating Imputation Paramters\n")
  condi <- mavis:::get_conditions(design)
  mis_vals <- data %>%
    dplyr::select(matches(condi)) %>%
    purrr::map(~which(is.na(.x)))

  mis_vals <- mis_vals[purrr::map_lgl(mis_vals, ~ length(.x) != 0)]

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
  dif <- vector()
  n_diff <- Inf
  row_imp <- function(row) {
    data <- row[-length(row)]
    if (!(all(is.na(data)) | !anyNA(data))) {
      data[is.na(data)] <- rnorm(
        sum(is.na(data)), mean(data, na.rm = T), row[length(row)]
      )
    }
    return(data)
  }

  imp_mat <- data %>%
    dplyr::select(matches(condi))

  tmp_na_vals <- imp_mat %>%
    dplyr::select(matches(condi)) %>%
    purrr::map(~which(is.na(.x)))

  tmp_na_vals <- tmp_na_vals[purrr::map_lgl(tmp_na_vals, ~ length(.x) != 0)]

  mean_vals <- imp_mat %>%
    purrr::map_dbl(mean, na.rm = T)
  for (i in names(tmp_na_vals)) {
    data[tmp_na_vals[[i]],i] <- imp_mat[tmp_na_vals[[i]],i] <- mean_vals[i]
  }

  while (TRUE) {
    tic <- Sys.time()
    for(i in names(mis_vals)) {
      form <- formula(
        paste0(i, '~ .')
      )
      idx <- mis_vals[[i]]

      pred <- ranger::ranger(
        form, dplyr::slice(imp_mat, -idx), mtry = floor(sum(design)/3),
        num.threads = workers
      ) %>%
        predict(dplyr::slice(imp_mat, idx)) %>%
        magrittr::use_series(predictions)

      dif[i] <- sum(
        (imp_mat[idx,i] - pred)^2
      ) / sum(pred^2)

      data[idx,i] <- imp_mat[idx,i] <- pred
    }
    if (n_diff > sum(dif)) {
      cat(
        'Previous error:', n_diff, '\t>\tCurrent error:', sum(dif), '\n'
      )

      n_diff <- sum(dif)
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
        trend_sd <- predict.glm(
          gam_reg,
          data.frame(mean = means[[i]]),
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

print_it_time <- function(tic) {
  toc <- Sys.time() - tic
  cat(
    "Iteration time:\n",
    format(toc), '\n\n'
  )
}
