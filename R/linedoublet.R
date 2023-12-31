
#' Create a line-doublet analytic element with a specified resistance
#'
#' [linedoublet()] creates a line-doublet analytic element with a constant specified resistance
#'
#' @param x0 numeric, starting x location of line-doublet
#' @param y0 numeric, starting y location of line-doublet
#' @param x1 numeric, ending x location of line-doublet
#' @param y1 numeric, ending y location of line-doublet
#' @param resistance numeric, hydraulic resistance of the line-doublet
#' @param ... ignored
#'
#' @details The resistance is used to calculate the constant strength of the line-doublet. If the aquifer is
#'    unconfined (i.e. has a variable saturated thickness), the system of equations becomes non-linear
#'    with respect to the hydraulic head and iteration is required to solve the model.
#'
#' For the special case that `resistance = Inf`, the line-doublet models an impermeable wall.
#'
#' @return Resistance-specified line-doublet analytic element which is an object of class `linedoublet` and inherits from `element`.
#' @export
#' @examples
#' linedoublet(-75, 50, 100, 50, res = 500)
#'
#' # Impermeable wall
#' linedoublet(-75, 50, 100, 50, resistance = Inf)
#'
linedoublet <- function(x0, y0, x1, y1, resistance, ...) {
  # TODO order, multiple control points
  order <- 0
  ld <- element(0, order + 1)
  ld$z0 <- x0 + y0 * 1i
  ld$z1 <- x1 + y1 * 1i
  ld$L <- abs(ld$z1 - ld$z0)

  ld$xc <- 0.5*(x0 + x1)
  ld$yc <- 0.5*(y0 + y1)
  ld$hc <- 0

  ld$order <- order
  if(resistance < 0) stop('Resistance should not be negative', call. = FALSE)
  ld$resistance <- resistance

  class(ld) <- c('linedoublet', class(ld))
  return(ld)
}

#'
#' @param linedoublet line-doublet analytic element of class `linedoublet` or inherits from it.
#'
#' @return complex potential influence of `linedoublet` evaluated at points `x y`.
#' @noRd
#'
omegainf.linedoublet <- function(linedoublet, x, y, ...) {
  zeta <- x + y * 1i
  Z <- (2 * zeta - (linedoublet$z0 + linedoublet$z1)) / (linedoublet$z1 - linedoublet$z0)
  tol <- 1e-12
  zp1 <- ifelse(abs(Z + 1) < tol, tol, Z + 1)
  zm1 <- ifelse(abs(Z - 1) < tol, tol, Z - 1)
  omi <- 1/(2*pi*1i) * log((zm1) / (zp1)) # order 0

  return(omi)
}

#'
#' @param linedoublet line-doublet analytic element of class `linedoublet` or inherits from it.
#'
#' @return complex discharge influence of `linedoublet` evaluated at points `x y`.
#' @noRd
#'
domegainf.linedoublet <- function(linedoublet, x, y, ...) {
  zeta <- x + y * 1i
  Z <- (2 * zeta - (linedoublet$z0 + linedoublet$z1)) / (linedoublet$z1 - linedoublet$z0)
  tol <- 1e-12
  zp1 <- ifelse(abs(Z + 1) < tol, tol, Z + 1)
  zm1 <- ifelse(abs(Z - 1) < tol, tol, Z - 1)
  m <- linedoublet$z1 - linedoublet$z0
  wi <- -(1/(2*pi*1i)) * zp1 * (2/(m*zp1) -  (2*zm1/(m*zp1^2))) / zm1 # order 0

  return(wi)
}

