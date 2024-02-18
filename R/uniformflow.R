
#' Create a uniform-flow analytic element with specified recharge
#'
#' [uniformflow()] creates a analytic element of constant uniform background flow.
#'
#' @param TR numeric, constant transmissivity value used to define the discharge.
#' @param gradient numeric, hydraulic gradient. Positive in the direction of flow.
#' @param angle numeric, angle of the primary direction of background flow,
#'    in degrees counterclockwise from the x-axis.
#' @param ... ignored
#'
#' @details `TR` and `gradient` are multiplied to obtain the discharge which remains
#'     constant throughout the system, independent of the saturated thickness of the aquifer.
#'
#' Groundwater flow is in the direction of the *negative* hydraulic gradient. Note that `gradient` is
#'     specified here as positive in the direction of flow for convenience.
#'
#' @return Analytic element of constant uniform flow which is an object of class `uniformflow` and inherits from `element`.
#' @export
#'
#' @examples
#' uniformflow(TR = 100, gradient = 0.002, angle = -45) # South-eastern direction
uniformflow <- function(TR, gradient, angle, ...) {
  uf <- element(gradient * TR)
  uf$udir <- exp(-1i * (angle*pi/180))
  class(uf) <- c('uniformflow', class(uf))
  return(uf)
}

#'
#' @param uf
#'
#' @return complex potential influence of `uf` evaluated at points `x y`.
#' @noRd
#'
omegainf.uniformflow <- function(uf, x, y, ...) {
  omi <- -uf$udir * (x + y * 1i)
  return(omi)
}

#'
#' @param uf
#'
#' @return complex discharge influence of `uf` evaluated at points `x y`.
#' @noRd
#'
domegainf.uniformflow <- function(uf, x, y, ...) {
  # handle NA when x = y = 0
  xx <- ifelse(x == 0 & y == 0, 1e-12, x)
  yy <- ifelse(x == 0 & y == 0, 1e-12, y)

  omi <- omegainf(uf, xx, yy, ...)
  zeta <- xx + yy*1i
  wi <- -omi/zeta
  return(wi)
}
