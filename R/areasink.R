
#' Create a circular area-sink analytic element with specified recharge
#'
#' [areasink()] creates a circular area-sink analytic element with constant, uniform specified recharge.
#'
#' @param xc numeric, x location of the center of the area-sink
#' @param yc numeric, y location of the center of the area-sink
#' @param N numeric, uniform constant leakage value (positive is into aquifer) in length per time.
#' @param R numeric, radius of the circular area-sink
#' @param location character, either `top` (default) or `base` specifying the vertical position of the area-sink.
#' @param ... ignored
#'
#' @details Area-sinks can be used to simulate areal recharge or seepage at the aquifer top, or leakage into
#'    or out of the aquifer at its base. The `location` argument is used when calculating the vertical flow
#'    component.
#'
#' @return Circular area-sink analytic element which is an object of class `areasink` and inherits from `element`.
#' @export
#'
#' @examples
#' areasink(xc = -500, yc = 0, N = 0.001, R = 500)
#'
#' dh <- 3
#' res <- 10 / 0.0001
#' areasink(xc = -500, yc = 0, N = -dh/res, R = 500, location = 'base')
#'
areasink <- function(xc, yc, N, R, location = c('top', 'base'), ...) {
  location <- match.arg(location)
  as <- element(N)
  as$xc <- xc
  as$yc <- yc
  as$R <- R
  as$location <- location
  class(as) <- c('areasink', class(as))
  return(as)
}

#'
#' @param areasink area-sink analytic element of class `areasink` or inherits from it.
#'
#' @return complex potential influence of `areasink` evaluated at points `x y`.
#' @noRd
#'
omegainf.areasink <- function(areasink, x, y, ...) {
  Rs <- areasink$R
  r <- sqrt((x - areasink$xc)^2 + (y - areasink$yc)^2)
  phi <- ifelse(r < Rs, -0.25*(r^2 - Rs^2), -Rs^2/2 * log(r/Rs)) + 0i
  return(phi)
}

#'
#' @param areasink area-sink analytic element of class `areasink` or inherits from it.
#'
#' @return complex discharge influence of `areasink` evaluated at points `x y`.
#' @noRd
#'
domegainf.areasink <- function(areasink, x, y, ...) {
  Rs <- areasink$R
  r <- sqrt((x - areasink$xc)^2 + (y - areasink$yc)^2)
  Qr <- ifelse(r <= Rs, r/2, Rs^2/(2*r)) # Haitjema 1995, eq. 5.30, for exfiltration
  W <- Qr*x/r - 1i*(Qr*y/r) # polar to cartesian
  W[r < 1e-15] <- 0 + 0i # continuous across boundary
  return(W)
}
