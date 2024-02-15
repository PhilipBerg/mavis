
# mavis

<!-- badges: start -->
[![R-CMD-check](https://github.com/PhilipBerg/mavis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilipBerg/mavis/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation
Note that `mavis` is currently only supported for `R` versions 4.2.X and to be compatible with RStan ([see here](https://blog.mc-stan.org/2022/04/26/stan-r-4-2-on-windows/)).
`mavis` requires `limma` to be installed; see [here](http://bioconductor.org/packages/release/bioc/html/limma.html).

You can install the development version of `mavis` as follows:

``` r
require(devtools)
devtools::install_github("PhilipBerg/mavis", build_vignettes = TRUE)
```

## Example

See the vignettes for examples:

``` r
library(mavis)
## See vignets for examples
vignette('baldur-tutorial') # For running Baldur
vignette('mult-imp') # For running multiple imputation and limma
```
