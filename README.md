
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mavis

<!-- badges: start -->

[![R-CMD-check](https://github.com/PhilipBerg/mavis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilipBerg/mavis/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

You can install the developmental version of `mavis` as follows:

``` r
require(devtools)
devtools::install_github("PhilipBerg/mavis", build_vignettes = TRUE)
```

## Example

See the vignettes for examples:

``` r
library(mavis)
## See vignettes for examples
vignette("Baldur-Tutorial")               # For running Baldur
vignette("Multiple Imputation and limma") # For running multiple imputation and limma
vignette("Ensemble Tutorial")             # For running the ensemble method
```
