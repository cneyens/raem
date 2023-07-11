
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raem

<!-- badges: start -->
<!-- badges: end -->

`raem` is an R package for modeling steady-state single-layer
groundwater flow using analytic elements.

## Installation

You can install the development version of raem from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cneyens/raem")
```

## Example

Construct an analytic element model with uniform background flow, two
extraction wells and a reference point.

Create elements:

``` r
library(raem)

TR = 100 # transmissivity
uf = uniformflow(TR, gradient = 0.001, angle = -45)
rf = constant(TR, xc = -1000, yc = 0, hc = 10)
w1 = well(xw = 200, yw = 0, Q = 250)
w2 = well(xw = -200, yw = -150, Q = 400)
```

Create the model. This automatically solves the system of equations.

``` r
m = aem(TR, uf, rf, w1, w2)
```

Plot head contours and element locations. First, specify the contouring
grid:

``` r
xg = seq(-500, 500, length = 100)
yg = seq(-250, 250, length = 100)
```

Now plot:

``` r
contour(m, xg, yg, z = 'head', col = 'dodgerblue3', nlevels = 20, asp = 1)
plot(m, add = TRUE)
```

<img src="man/figures/README-plot-head-1.png" width="80%" />

Compute particle traces along `y = 200` at 20 intervals per year for 10
years and add to plot:

``` r
paths = tracelines(m, x0 = seq(-450, 450, 50), y0 = 200, times = seq(0, 10*365, 365/20))
plot(paths, add = TRUE, col = 'orange3')
```

<img src="man/figures/README-plot-traces-1.png" width="80%" />
