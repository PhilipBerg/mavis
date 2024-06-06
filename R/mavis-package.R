#' @keywords internal
"_PACKAGE"

#' The 'mavis' package.
#'
#' @description
#' Mavis is completely data-driven from the observed quantifications and can be applied to any post-quantification data set.
#' Mavis can perform both imputation and statistical decisions.
#' For imputation, Mavis produces a more robust estimator of the variance component of the missing data by using the M-V trend.
#' Then, for the statistical testing for mean differences, Mavis makes inferences on the individual data points using the M-V trend.
#' This, can either be used with the Baldur framework or the limma framework.
#'
#' @name mavis-package
#' @aliases mavis
#' @import methods
#' @import Rcpp
## usethis namespace: start
#' @importFrom baldur empirical_bayes
#' @importFrom baldur fit_gamma_regression
#' @importFrom baldur plot_gamma
#' @importFrom baldur psrn
#' @importFrom baldur weakly_informative
## usethis namespace: end
NULL
