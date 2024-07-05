
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raem

<!-- badges: start -->

[![R-CMD-check](https://github.com/cneyens/raem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cneyens/raem/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`raem` is an R package for modeling steady-state single-layer
groundwater flow under the Dupuit-Forchheimer assumption using analytic
elements.

## Installation

You can install the development version of `raem` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cneyens/raem")
```

## Example

Construct an analytic element model of an aquifer with uniform
background flow, two extraction wells and a reference point.

Specify the aquifer parameters and create elements:

``` r
library(raem)

k = 10
top = 10
base = 0
n = 0.2

hr = 15 # head at reference point, above aquifer top
TR = k * (top - base) # constant transmissivity of background flow

uf = uniformflow(TR, gradient = 0.001, angle = -45)
rf = constant(xc = -1000, yc = 0, hc = hr)
w1 = well(xw = 200, yw = 0, Q = 250)
w2 = well(xw = -200, yw = -150, Q = 400)
```

Create the model. This automatically solves the system of equations.

``` r
m = aem(k = k, top = top, base = base, n = n, 
        uf, rf, w1, w2)
```

Find the head and discharge at two locations `x = -200, y = 200` and
`x = 100, y = 200`. Note that there are no vertical flow components in
this model:

``` r
heads(m, x = c(-200, 100), y = 200)
#> [1] 13.64573 13.33314
```

``` r

discharge(m, c(-200, 100), 200, z = top) # [L^2 / T]
#>              Qx         Qy Qz
#> [1,] 0.15028815 -0.2923908  0
#> [2,] 0.06041242 -0.3347206  0
```

Plot head contours and element locations. First, specify the contouring
grid:

``` r
xg = seq(-500, 500, length = 100)
yg = seq(-250, 250, length = 100)
```

Now plot:

``` r
contours(m, xg, yg, 'heads', col = 'dodgerblue', nlevels = 20)
plot(m, add = TRUE)
```

<img src="man/figures/README-plot-head-1.png" width="100%" />

Compute particle traces starting along `y = 200` at 20 intervals per
year for 5 years and add to plot:

``` r
paths = tracelines(m, x0 = seq(-450, 450, 50), y0 = 200, z0 = top, times = seq(0, 5*365, 365/20))
plot(paths, add = TRUE, col = 'orange')
```

<img src="man/figures/README-plot-traces-1.png" width="100%" />
