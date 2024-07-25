
#' Create a constant-head analytic element
#'
#' [constant()] creates an analytic element containing a constant head, often referred to
#'     as *the reference point*.
#'
#' @param xc numeric, x location of the reference point.
#' @param yc numeric, y location of the reference point.
#' @param hc numeric, hydraulic head at the reference point.
#' @param ... ignored
#'
#' @return Constant-head analytic element point which is an object of class `constant` and inherits from `element`.
#' @export
#'
#' @examples
#' rf <- constant(xc = -100, yc = 0, hc = 10)
#'
constant <- function(xc, yc, hc, ...) {
  cn <- element(0, 1)
  cn$xc <- xc
  cn$yc <- yc
  cn$hc <- hc
  class(cn) <- c('constant', class(cn))
  return(cn)
}

#'
#' @param constant an constant-head analytic element of class `constant`.
#'
#' @return complex potential influence of `constant` evaluated at points `x y`.
#' @noRd
#'
omegainf.constant <- function(constant, x, y, ...) {
  crds <- cbind(x, y)
  omi <- complex(length.out = nrow(crds), real = 1, imaginary = 0)
  return(omi)
}

#'
#' @param constant an constant-head analytic element of class `constant`.
#'
#' @return complex discharge influence of `constant` evaluated at points `x y`.
#' @noRd
#'
domegainf.constant <- function(constant, x, y, ...) {
  crds <- cbind(x, y)
  wi <- complex(length.out = nrow(crds), real = 0, imaginary = 0)
  return(wi)
}
