
#' Create a uniform-flow analytic element with specified recharge
#'
#' [uniformflow()] creates a analytic element of constant uniform background flow.
#'
#' @param TR numeric, transmissivity value of aquifer
#' @param gradient numeric, hydraulic gradient. Positive in the direction of flow.
#' @param angle numeric, angle of the primary direction of background flow,
#'    in degrees counterclockwise from the x-axis.
#' @param ... ignored
#'
#' @return Analytic element of constant uniform flow which is an object of class `uf` and inherits from `element`.
#' @export
#'
#' @examples
#' uniformflow(gradient = 0.002, angle = -45, TR = 100) # South-eastern direction
uniformflow <- function(TR, gradient, angle, ...) {
  uf <- element(gradient * TR)
  uf$udir <- exp(-1i * (angle*pi/180))
  class(uf) <- c('uf', class(uf))
  return(uf)
}

#'
#' @param uf
#'
#' @return complex potential influence of `uf` evaluated at points `x y`.
#' @noRd
#'
omegainf.uf <- function(uf, x, y, ...) {
  omi <- -uf$udir * (x + y * 1i)
  return(omi)
}

#'
#' @param uf
#'
#' @return complex discharge influence of `uf` evaluated at points `x y`.
#' @noRd
#'
domegainf.uf <- function(uf, x, y, ...) {
  omi <- omegainf(uf, x, y, ...)
  zeta <- x + y*1i
  wi <- -omi/zeta
  return(wi)
}
