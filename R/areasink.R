
#' Create a circular area-sink analytic element with specified recharge
#'
#' [areasink()] creates a circular area-sink analytic element with constant, uniform specified recharge.
#'
#' @param xc numeric, x location of the center of the area-sink
#' @param yc numeric, y location of the center of the area-sink
#' @param N numeric, uniform constant recharge value (positive is into aquifer) in length per time.
#' @param R numeric, radius of the circular area-sink
#' @param ... ignored
#'
#' @return Circular area-sink analytic element which is an object of class `areasink` and inherits from `element`.
#' @export
#'
#' @examples
#' areasink(xc = -500, yc = 0, N = 0.001, R = 500)
#'
areasink <- function(xc, yc, N, R, ...) {
  as <- element(N)
  as$xc <- xc
  as$yc <- yc
  as$R <- R
  class(as) <- c('areasink', class(as))
  return(as)
}

#'
#' @param areasink
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

# TODO fix this
#'
#' @param areasink
#'
#' @return complex discharge influence of `areasink` evaluated at points `x y`.
#' @rdname domega
#'
domega.areasink <- function(areasink, x, y, ...) { # TODO, absolute values not correct
  # Haitjema 1995, eq. 5.30, for exfiltration
  Rs <- areasink$R
  r <- sqrt((x - areasink$xc)^2 + (y - areasink$yc)^2)
  # Bakker & Post, 2022, eq. 6.39 & 6.40
  Qr <- areasink$parameter * ifelse(r <= Rs, r/2, Rs^2/(2*r))
  W <- (Qr*(x - areasink$xc)/(2*pi*r^2)) + (-Qr*(y - areasink$yc)/(2*pi*r^2))*1i # not sure about this
  W[r < 1e-15] <- 0 + 0i # continuous across boundary
  return(W)
}
