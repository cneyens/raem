
#' Create a analytic element of a constant-discharge well
#'
#' [well()] creates an analytic element of a well with constant discharge.
#'
#' @param xw numeric, x location of the well
#' @param yw numeric, y location of the well
#' @param Q numeric, volumetric discharge of the well (positive is out of aquifer)
#' @param rw numeric, radius of well. Defaults to 0.3 (meter).
#' @param ... ignored
#'
#' @return Analytic element of a well with constant discharge which is an object of class `well` and inherits from `element`.
#' @export
#'
#' @examples
#' well(xw = 50, yw = 0, Q = 200, rw = 0.3)
#'
well <- function(xw, yw, Q, rw = 0.3, ...) {
  well <- element(Q)
  well$zetaw <- xw + 1i * yw
  well$rw <- rw
  class(well) <- c('well', class(well))
  return(well)
}

#' Create a analytic element of a well with a constant head
#'
#' [headwell()] creates an analytic element of a well with a constant, specified head. The discharge
#'    into the well is computed by solving the corresponding `aem` model.
#'
#' @param xw numeric, x location of the well
#' @param yw numeric, y location of the well
#' @param hw numeric, specified hydraulic head in the well
#' @param rw numeric, radius of the well. Defaults to 0.3 (meter).
#' @param ... ignored
#'
#' @return Analytic element of a well with constant head which is an object of class `headwell` and inherits from `well`.
#' @export
#'
#' @examples
#' headwell(xw = 400, yw = 300, hw = 20, rw = 0.3)
#'
headwell <- function(xw, yw, hw, rw = 0.3, ...) {
  hwe <- well(xw = xw, yw = yw, Q = 0, rw = rw)
  hwe$xc <- xw + rw
  hwe$yc <- yw
  hwe$hc <- hw
  hwe$nunknowns <- 1
  class(hwe) <- c('headwell', 'headequation', class(hwe))
  return(hwe)
}

#'
#' @param well well analytic element of class `well` or inherits from it.
#'
#' @return complex potential influence of `well` evaluated at points `x y`.
#' @noRd
#'
omegainf.well <- function(well, x, y, ...) {
  zminzw <- (x + 1i*y) - well$zetaw
  alpha <- atan2(y, x)
  zminzw <- ifelse(abs(zminzw) < well$rw, well$rw*exp(alpha*1i), zminzw)
  omi <- 1/(2*pi) * log(zminzw)
  return(omi)
}

#'
#' @param well well analytic element of class `well` or inherits from it.
#'
#' @return complex discharge influence of `well` evaluated at points `x y`.
#' @noRd
#'
domegainf.well <- function(well, x, y, ...) {
  zminzw <- (x + 1i*y) - well$zetaw
  alpha <- atan2(y, x)
  zminzw <- ifelse(abs(zminzw) < well$rw, well$rw*exp(alpha*1i), zminzw)
  wi <- -1/(2*pi*zminzw)
  return(wi)
}
