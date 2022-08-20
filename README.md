
# mavis

<!-- badges: start -->
<!-- badges: end -->

## Installation
Note that `mavis` is currently only supported for `R` versions 4.0.X and will be updated once RStan is compatible ([see here](https://blog.mc-stan.org/2022/04/26/stan-r-4-2-on-windows/)).
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
