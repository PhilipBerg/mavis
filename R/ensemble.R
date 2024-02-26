#' @title An R6 class for Bayesian Ensemble
#'
#' @description
#' An [R6::R6Class] to interface with the Bayesian ensemble.
#'
#' @importFrom R6 R6Class
#' @importFrom utils find
#' @docType class
#'
#' @field stack A data stack for the methods to have available for the use in
#'    the ensemble.
#' @name ensemble
#' @inherit start_ensemble examples
ensemble <- R6::R6Class("Ensemble_Stack",
                        list(
                          stack = list(),
                          #' @param data Input data with the results from a
                          #'    method.
                          #' @param method The name to use to refer to the
                          #'    method in the stack.
                          #' @param id_col,p_col,lfc_col,comp_col  Name of the
                          #'    columns containing the unique identifier for
                          #'    each peptide, the p-values, log-fold change, and
                          #'    contast/comparison of the method, respectively.
                          #' @param auxilary Additional paramters to store in
                          #'    the stack. Can be useful if `do_rm` is `TRUE`.
                          #'    Defaults to "rest" meaning the rest of the
                          #'    columns (other than the forementioned).
                          #' @param do_rm Should the variable be removed from
                          #'    its current environment? Can be useful to reduce
                          #'    memory usage and avoid the need to rm after each
                          #'    method is added to the stack. `ensemble` will
                          #'    first try to find the environment of the `data`,
                          #'    if unsuccessful, it will try to remove it from
                          #'    the caller environment. Defaults to `FALSE`.
                          #'
                          #' @description
                          #' Add a new method to the stack
                          add   = function(data, method, id_col, p_col, lfc_col, comp_col = "comparison", auxilary = "rest", do_rm = FALSE) {
                            dq <- substitute(data)
                            self$stack[[method]]$id_col   <- id_col
                            self$stack[[method]]$p_col    <- p_col
                            self$stack[[method]]$lfc_col  <- lfc_col
                            self$stack[[method]]$comp_col <- comp_col
                            if (auxilary == "rest") {
                              self$stack[[method]]$data <- data
                            } else if (auxilary == "none") {
                              self$stack[[method]]$data <- data[c(id_col, p_col, lfc_col, comp_col)]
                            } else {
                              self$stack[[method]]$data <- data[c(id_col, p_col, lfc_col, comp_col, auxilary)]
                            }
                            if (do_rm) {
                              ev <- find(as.character(dq))
                              if (rlang::is_empty(ev)) {
                                ev <- rlang::caller_env()
                              } else {
                                ev  <- get(ev)
                              }
                              do.call("rm", list(dq, envir = ev))
                            }
                            invisible(self)
                          },
                          #' @description
                          #' Remove a method from the stack and return it to stdout
                          #'
                          #' @param method method index (integer) or name to return
                          #' defaults to the most recently added method.
                          pop  = function(method = length(self$stack)) {
                            if (length(method) == 1) {
                               f    <- `[[`
                              `f<-` <- `[[<-`
                            } else {
                               f    <- `[`
                              `f<-` <- `[<-`
                            }
                            out <- f(self$stack, method)
                            f(self$stack, method) <- NULL
                            out
                          },
                          #' @description
                          #' Run the ensemble
                          #'
                          #' @param methods Which methods to include in the
                          #'    ensemble.
                          #' @param parallel_chains Number of parallel workers
                          #'    to use when running the sampling.
                          #' @param parallel_runs The number of parallel
                          #'    comparison to run at the same time. Note that
                          #'    the total number of parallel instances will be
                          #'    `parallel_chains`\eqn{\times}`parallel_runs`.
                          #'    Can be efficient when there are many comparisons
                          #'    to run.
                          #' @param ... Additional parameters to
                          #'    [rstan::sampling].
                          #'
                          #' @return A [tibble::tibble] with the ensembled
                          #'    p-value. Column names will default to the first
                          #'    method in the stack.
                          run_ensemble = function(methods = "all", parallel_chains = 1, parallel_runs = 1, ...) {
                            if (length(methods) == 1) if(methods == "all") methods <- names(self$stack)
                            col_names <- self$stack[methods][[1]][1:4]
                            stn_data <- private$make_stan_inputs(methods)
                            id_order <- private$standardize_data(methods, simplify = TRUE)$id
                            stan_mod <- stanmodels$mavis

                            stn_inputs <- rlang::dots_list(...)
                            stn_call <- rlang::expr(
                              ~ rstan::sampling(
                                object = stan_mod,
                                cores = parallel_chains,
                                data = .x,
                                !!!stn_inputs
                              )
                            )
                            if (parallel_runs == 1) {
                              samp <- purrr::map(stn_data,
                                                 rlang::eval_tidy(stn_call)
                              )
                            } else {
                              cl <- multidplyr::new_cluster(parallel_runs)
                              suppressWarnings(invisible((
                                purrr::quietly(multidplyr::cluster_library))(cl,
                                                                             c(
                                                                               "dplyr", "tidyr", "purrr", "tibble", "stringr",
                                                                               "magrittr", "rstan", "StanHeaders", "rlang"
                                                                             )
                                ))
                              )
                              multidplyr::cluster_copy(cl, "stn_call")
                              multidplyr::cluster_copy(cl, "stan_mod")
                              samp <- tibble::enframe(stn_data) %>%
                                multidplyr::partition(cl) %>%
                                dplyr::mutate(
                                  samp = purrr::map(value, rlang::eval_tidy(stn_call))
                                ) %>%
                                dplyr::collect()
                            }
                            samp %>%
                              private$stan_to_output() %>%
                              dplyr::mutate(
                                !!col_names$id_col := id_order,
                                mean = 1 - mean
                              ) %>%
                              dplyr::select(all_of(col_names$id_col), !!col_names$comp_col := name, p_val = mean, everything())
                          },
                          #' @description
                          #' Printing method for the class.
                          print = function() {
                            n <- length(self$stack)
                            top <- cli::cli_text("Ensemble Stack of {n}")
                            if(n > 0) {
                              cli::cli_text("Methods Available: {paste0(names(self$stack), collapse = ', ')}")
                            }
                            invisible(self)
                          }
                        )
)

#' Start a New Bayesian Ensemble
#'
#' @description
#' Creates a new ensemble instance.
#'
#' @name start_ensemble
#'
#' @return An ensemble object; see [ensemble] for details.
#'
#' @examples
#' # Setup variables
# Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#' contrast <- limma::makeContrasts(
#'   contrasts = "ng100-ng50",
#'   levels = design
#' )
#'
#' # Exemplify on the non-missing data
#' yeast <- tidyr::drop_na(yeast) %>%
#'    # Normalize and log-transform the data
#'    psrn("identifier")
#'
#' # Fit the gamma regressions
#' gamma_reg_model <- yeast %>%
#'    calculate_mean_sd_trends(design) %>%
#'    fit_gamma_regression(formula = sd ~ mean)
#'
#' limma_gr <- run_limma(
#'   yeast,
#'   design,
#'   contrast,
#'   gamma_reg_model,
#'   "identifier"
#' ) %>%
#'    head(10)
#' limma_trend <- run_limma(
#'   yeast,
#'   design,
#'   contrast,
#'   NULL,
#'   "identifier"
#' ) %>%
#'    dplyr::filter(identifier %in% limma_gr$identifier)
#'
#' new_ensemble <- start_ensemble()
#' new_ensemble$add(limma_gr, "GR-limma", "identifier", "p_val", "lfc", do_rm = TRUE
#' )$add(
#'    limma_trend, "limma-trend", "identifier", "p_val", "lfc", do_rm = TRUE
#' )
#' ensemble_results <- new_ensemble$run_ensemble(iter = 150)
#' limma_gr <- new_ensemble$pop("limma_gr")
#' @export
start_ensemble <- function() {
  ensemble$new()
}

ensemble$set(
  "private", "std_data",
  tibble::tibble()
)

ensemble$set(
  "private", "stan_to_output",
  function(sample) {
    if (tibble::is_tibble(sample)) {
      sample <- dplyr::select(sample, -value)
    } else {
      sample <- tibble::enframe(sample, value = "samp")
    }
    sample %>%
      dplyr::mutate(
        samp = purrr::map(samp, ~ .x %>%
                            rstan::summary(pars = "mu") %>%
                            magrittr::use_series(summary) %>%
                            tibble::as_tibble(rownames = 'par') %>%
                            dplyr::select(mean, n_eff, Rhat)
        )
      ) %>%
      tidyr::unnest(cols = samp)
  }
)

ensemble$set(
  "private", "standardize_data",
  function(methods, simplify) {
    out <- self$stack[methods] %>%
      purrr::imap(
        ~ dplyr::select(.x$data,
                        id = !!.x$id_col,
                        !!paste0(.y, "_pv") := !!.x$p_col,
                        !!paste0(.y, "_lfc") := !!.x$lfc_col,
                        comp = comparison
        )
      )
    if (simplify) {
      purrr::reduce(out, dplyr::left_join, by = dplyr::join_by(id, comp))
    } else {
      out
    }
  }
)

ensemble$set(
  "private", "make_stan_inputs", function(methods) {
    std_data <- private$standardize_data(methods, simplify = FALSE)
    weights  <- private$make_weights(std_data)
    purrr::reduce(std_data, dplyr::left_join, by = dplyr::join_by(id, comp)) %>%
      split.data.frame(.$comp) %>%
      purrr::map(
        private$make_stan_input
      ) %>%
      purrr::map2(weights, private$add_weights)
  }
)

ensemble$set(
  "private", "make_stan_input", function(decision) {
    N <- nrow(decision)
    pvals <- decision %>%
      dplyr::select(matches("_pv$")) %>%
      dplyr::mutate(
        across(where(is.numeric), ~ dplyr::if_else(. != 0, ., min(.[.!=0])))
      )
    M <- ncol(pvals)
    list(
      N = N,
      M = M,
      y =
        matrix(
          unlist(pvals), N, M
        )
    )
  }
)

ensemble$set(
  "private", "make_weights", function(std_data) {
    weights <- std_data %>%
      purrr::map(~ split.data.frame(.x, .x$comp)) %>%
      purrr::map_depth(2,
                       ~ cor(.x[[2]], abs(.x[[3]]), method = "spear")
      ) %>%
      purrr::map(unlist)

    weights <- names(weights[[1]]) %>%
      setNames(., .) %>%
      purrr::map(
        ~ lapply(weights, `[[`, .x)
      ) %>%
      purrr::map(unlist) %>%
      purrr::map(
        ~ .x/sum(.x)
      )
    weights
  }
)

ensemble$set(
  "private", "add_weights", function(stan_input, weights) {
    stan_input$tau <- weights
    stan_input
  }
)
