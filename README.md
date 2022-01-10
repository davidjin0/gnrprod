
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gnrprod`

`gnrprod` implements the nonparametric identification of gross output
production functions outlined in Gandhi, Navarro, and Rivers (2020). The
main wrapper function `gnrprod` estimates the parameters of production
functions and productivity. The current version of `gnrprod` supports
only one flexible input.

## Installation

To install the latest development version of  from Github, run the
following from R:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("davidjin0/gnrprod")
library(gnrprod)
```

## Usage

This package features three main functions: `gnrprod`, `gnrflex`, and
`gnriv`. `gnrprod` performs the entire production function estimation
routine. `gnrflex` performs the first stage of the estimation routine
which returns the flexible input elasticity. `gnrflex` performs the
second stage of the estimation routine which returns the fixed input
elasticities and productivity.

An example of the use of `gnrprod`:

``` r
require(gnrprod)
#> Loading required package: gnrprod

# load Colombian plant-level data
data <- colombian

# estimate production function parameters and productivity
gnr_fit <- gnrprod(output = "RGO", fixed = c("L", "K"), flex = "RI",
                   share = "share", id = "id", time = "year", data = data)

# print results
# gnr_fit
# summary(gnr_fit)
```

Alternatively, one can use `gnrflex` and `gnriv`:

``` r
# estimate flexible input elasticities (first stage)
gnr_fs <- gnrflex(output = "RGO", fixed = c("L", "K"), flex = "RI",
                  share = "share", id = "id", time = "year", data = data)

# print estimate
# gnr_fs

# estimate second stage
gnr_ss <- gnriv(object = gnr_fs)

# print estimates
# gnr_ss
```

## References

Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. “On the
Identification of Gross Output Production Functions.” *Journal of
Political Economy*, 128(8): 2973-3016. <https://doi.org/10.1086/707736>.
