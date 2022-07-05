#' The 'mavis' package.
#'
#' @description
#' Mavis is completely data-driven from the observed quantifications and can be applied to any post-quantification data set.
#' Mavis can perform both imputation and statistical decisions.
#' For imputation, Mavis produces a more robust estimator of the variance component of the missing data by using the M-V trend.
#' Then, for the statistical testing for mean differences, Mavis makes inferences on the individual data points using the M-V trend.
#' This, can either be used with the Baldur framework or the limma framework.
#'
#' @docType package
#' @name mavis-package
#' @aliases mavis
#' @useDynLib mavis, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.5. https://mc-stan.org
#'
NULL