
#' Create a analytic element of a constant-discharge well
#'
#' [well()] creates an analytic element of a well with constant discharge.
#'
#' @param xw numeric, x location of the well
#' @param yw numeric, y location of the well
#' @param Q numeric, volumetric discharge of the well (positive is out of aquifer)
#' @param rw numeric, radius of well. Defaults to 0.3.
#' @param ... ignored
#'
#' @return Analytic element of a well with constant discharge which is an object of class `well` and inherits from `element`.
#' @export
#' @seealso [headwell()]
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
#'    into the well is computed by solving the corresponding `aem` model. The head can be specified at
#'    the well or at any other location.
#'
#' @param xw numeric, x location of the well
#' @param yw numeric, y location of the well
#' @param hc numeric, specified hydraulic head at the collocation point.
#' @param rw numeric, radius of the well. Defaults to 0.3 (meter).
#' @param xc numeric, x location of the collocation point. See details. Defaults to `xw`.
#' @param yc numeric, y location of the collocation point. See details. Defaults to `yw`.
#' @param rc numeric, radius of the collocation point. See details. Defaults to `rw`.
#' @param resistance numeric, hydraulic resistance of the well screen at `xw yw`. Defaults to 0 (no resistance).
#' @param ... ignored
#'
#' @details The discharge from the well at location `xw yw` is computed by solving the `aem` model given
#'    the specified head `hc`. This head can be specified at any point, called the collocation point.
#'    This can be used to compute the discharge of the well by specifying the head at some other location.
#'    The head is specified at `xc + rc`, `yc`. By default, the location of the well and the collocation point are the same.
#'
#' The resistance at the well screen can be increased for a well in poor connection with the aquifer due. If the aquifer is
#'    unconfined (i.e. has a variable saturated thickness), the system of equations becomes non-linear
#'    with respect to the hydraulic head and iteration is required to solve the model.
#'
#' @return Analytic element of a well with constant head which is an object of class `headwell` and inherits from `well`.
#' @export
#' @seealso [well()]
#' @examples
#' headwell(xw = 400, yw = 300, hc = 20, rw = 0.3)
#' headwell(xw = 400, yw = 300, hc = 20, rw = 0.3, resistance = 10)
#' headwell(xw = 400, yw = 300, hc = 20, rw = 0.3, xc = 500, yc = 500, rc = 0)
#'
headwell <- function(xw, yw, hc, rw = 0.3, xc = xw, yc = yw, rc = rw, resistance = 0, ...) {
  hwe <- well(xw = xw, yw = yw, Q = 0, rw = rw)
  hwe$xc <- xc + rc
  hwe$yc <- yc
  hwe$hc <- hc
  hwe$nunknowns <- 1
  hwe$resistance <- resistance
  class(hwe) <- c('headwell', class(hwe))
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
