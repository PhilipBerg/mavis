#' Add two variables together
#'
#' This function adds the values of x and y.
#'
#' @param x Numeric value for one of the variables
#' @param y Numeric value of the other variable
#'
#' @name add
#' @return A numeric value of the sum of x and y
#' @export
#'
#' @examples
#' # Common part
#' x <- 1
#' y <- 2
#' # Unique part compared to product
#' z <- add(x, y)
add <- function(x, y) {
  x + y
}

#' Multiply two variables together
#'
#' This function takes the product of x and y
#' @inheritParams add
#'
#' @return A numeric value of the product of x and y
#' @export
#'
#' @examples
#' knitr::knit_child("man/rmd/child_test.Rmd")
#'
#' # Unique part
#' product(x, y)
product <- function(x, y) {
  x * y
}


#' Calculates the mean of the summed variable x
#'
#' @param x Sum of variables to calculate the mean from
#' @param n Number of variables that were summed
#'
#' @return A numerical value for the mean of x
#' @export
#'
#' @examples
#' # Common
#' x <- 1
#' y <- 2
#' z <- add(x, y)
#' # Unique
#' n <- 2
#' mean2(z, n)
mean2 <- function(x, n) {
  x*n^-1
}
