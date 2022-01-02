
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RPSHA

<!-- badges: start -->

[![GitHub
license](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/wltcwpf/RPSHA/blob/master/LICENSE.md)
[![Report
Issues!](https://img.shields.io/badge/Report%20Issues-Here-1abc9c.svg)](https://github.com/wltcwpf/RPSHA/issues)
[![Open Source?
Yes!](https://img.shields.io/badge/Open%20Source-Yes-green.svg)](https://github.com/wltcwpf/RPSHA)
<!-- badges: end -->

RPSHA stands for Regional Probabilistic Seismic Hazard Assessment. This
package is developed to conduct the conventional PSHA for independent
sites in California using UCERF3 source model. More importantly, this
package implements LASSO regression for event selection and ground
motion map selection to assess seismic hazard for regional distributed
infrastructures that are no longer independent but spatially
distributed.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
devtools::install_github("wltcwpf/RPSHA")
```
