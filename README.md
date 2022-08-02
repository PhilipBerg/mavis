
# mavis

<!-- badges: start -->
<!-- badges: end -->

## Installation
Note that `mavis` is currently only supported for `R` versions 4.0.X 
You can install the development version of `mavis` as follows:

``` r
require(devtools)
devtools::install_github("PhilipBerg/mavis", build_vignettes = TRUE)
```

## Example

See the vignets for examples:

``` r
library(mavis)
## See vignets for examples
vignette('baldur-tutorial') # For running Baldur
vignette('mult-imp') # For running multiple imputation and limma
```
