
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
#' @details The resistance is used to calculate the strength of the line-doublet at the collocation points.
#'  If the aquifer is unconfined (i.e. has a variable saturated thickness), the system of equations becomes
#'  non-linear with respect to the hydraulic head and iteration is required to solve the model.
#'
#' For the special case that `resistance = Inf`, the line-doublet models an impermeable wall.
#'
#' Three collocation points are defined per line-doublet: one at each vertex of the line-doublet and one in
#'  the center of the line. A single line-doublet therefore adds three equations to the system of equations.
#'  The strength distribution is defined at these collocation points and is approximated by a parabolically varying
#'  strength distribution in between.
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
  order <- 2 # only parabolic strength allowed
  ld <- element(rep(0, order + 1), order + 1)
  ld$z0 <- x0 + y0 * 1i
  ld$z1 <- x1 + y1 * 1i
  ld$L <- abs(ld$z1 - ld$z0)

  ld$xc <- c(x0, 0.5*(x1 - x0) + x0, x1)
  ld$yc <- c(y0, 0.5*(y1 - y0) + y0, y1)

  # Chebyshev roots # TODO tests fail at collocation points except at center
  # cheby_root <- cos(pi * (seq(order-1) + 0.5) / order)
  # half_length <- 0.5*(ld$z1 - ld$z0)
  # ld$xc <- c(cheby_root * Re(half_length), 0.5*(x0 + x1), -cheby_root * Re(half_length))
  # ld$yc <- c(cheby_root * Im(half_length), 0.5*(y0 + y1), -cheby_root * Im(half_length))

  ld$hc <- rep(0, order + 1)

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

  # parabolic strength
  fm <- -0.5 * zm1 * log(zm1 / zp1) - 1
  gm <- 0.5 * zp1 * log(zm1 / zp1) + 1
  pm <- (Z^2 - 1) * log(zm1 / zp1) + 2*Z

  # return in matrix form so omega.element can perform a matrix multiplication with the parameter vector
  omi <- 1/(2*pi*1i) * cbind(fm + pm/2, -pm, gm + pm/2)

  # omi <- 1/(2*pi*1i) * (s1 * fm + s2*gm + s1*pm/2 + s2*pm/2 - sc*pm) # parabolic, order=2
  # omi <- 1/(2*pi*1i) * log((zm1) / (zp1)) # constant strength, order=0; deprecated

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
  Znum <- (2 * zeta - (linedoublet$z0 + linedoublet$z1))

  # parabolic strength
  fm <- -log(zm1 / zp1) / m - 0.5 * (zp1 * (2 / (m * zp1) - (2*zm1) / (m * zp1^2)))
  gm <- log(zm1 / zp1) / m + (zp1^2) * (2 / (m * zp1) - 2 * zm1 / (m * zp1^2)) / (2 * zm1)
  pm <- ((4 * Znum * log(zm1 / zp1)) / m^2) + (zp1 * (Z^2 - 1) * (2 / (m * zp1) - 2 * zm1 / (m * zp1^2))) / zm1

  # return in matrix form so omega.element can perform a matrix multiplication with the parameter vector
  wi <- -1/(2 * pi * 1i) * cbind(fm + pm/2, -pm, gm + pm/2)

  # wi <- -(1/(2*pi*1i)) * zp1 * (2/(m*zp1) -  (2*zm1/(m*zp1^2))) / zm1 # order 0; deprecated

  return(wi)
}

