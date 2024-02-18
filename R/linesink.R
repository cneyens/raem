
#' Create a strength-specified line-sink analytic element
#'
#' [linesink()] creates a line-sink analytic element with constant specified strength.
#'
#' @param x0 numeric, starting x location of line-sink
#' @param y0 numeric, starting y location of line-sink
#' @param x1 numeric, ending x location of line-sink
#' @param y1 numeric, ending y location of line-sink
#' @param sigma numeric, specific strength of the line-sink, i.e. discharge per length of line-sink. Positive is out of aquifer.
#' @param width numeric, width of the line-sink. Only used in [tracelines()] to determine if a particle has reached the line. Defaults to zero (infinitesimally narrow line).
#' @param ... ignored
#'
#' @return Strength-specified line-sink analytic element which is an object of class `linesink` and inherits from `element`.
#' @export
#' @seealso [headlinesink()]
#' @examples
#' linesink(-75, 50, 100, 50, sigma = 1)
#'
linesink <- function(x0, y0, x1, y1, sigma, width = 0, ...) {
  ls <- element(sigma)
  ls$z0 <- x0 + y0 * 1i
  ls$z1 <- x1 + y1 * 1i
  ls$L <- abs(ls$z1 - ls$z0)
  ls$width <- width
  class(ls) <- c('linesink', class(ls))
  return(ls)
}

#' Create a head-specified line-sink analytic element
#'
#' [headlinesink()] creates a line-sink analytic element with constant specified head. The discharge
#'    into the line-sink per unit length is computed by solving the corresponding `aem` model.
#'
#' @param x0 numeric, starting x location of line-sink
#' @param y0 numeric, starting y location of line-sink
#' @param x1 numeric, ending x location of line-sink
#' @param y1 numeric, ending y location of line-sink
#' @param hc numeric, specified hydraulic head of the line-sink
#' @param resistance numeric, hydraulic resistance of the line-sink at its connection with the aquifer. Defaults to 0 (no resistance).
#' @param width numeric, width of the line-sink. Used with `resistance` to calculate the line-sink strength, and by [tracelines()] to
#'    determine if a particle has reached the line. Defaults to zero (infinitesimally narrow line).
#' @param ... ignored
#'
#' @details The strength of the line-sink (discharge per unit length of line-sink) is computed by solving
#'    the `aem` model given the specified head `hc` for the line-sink. This head is located at the so-called
#'    collocation point, which is placed at the center of the line-sink.
#'
#' The resistance can be increased for a line-sink in poor connection with the aquifer. The effect of a larger
#'    or smaller wetted perimeter can be mimicked by adjusting the `resistance` and/or `width` accordingly. If
#'    `width = 0` (the default) it is removed from the conductance calculation. If the aquifer is unconfined
#'    (i.e. has a variable saturated thickness), the system of equations becomes non-linear with respect to the
#'    hydraulic head and iteration is required to solve the model.
#'
#' @return Head-specified line-sink analytic element which is an object of class `headlinesink` and inherits from `linesink`.
#' @export
#' @seealso [linesink()]
#' @examples
#' headlinesink(-75, 50, 100, 50, hc = 10)
#' headlinesink(-75, 50, 100, 50, hc = 10, res = 10, width = 4)
#'
headlinesink <- function(x0, y0, x1, y1, hc, resistance = 0, width = 0, ...) {
  hls <- linesink(x0 = x0, y0 = y0, x1 = x1, y1 = y1, sigma = 0, width = width)
  hls$xc <- 0.5*(x0 + x1)
  hls$yc <- 0.5*(y0 + y1)
  hls$hc <- hc
  hls$nunknowns <- 1
  hls$resistance <- resistance
  class(hls) <- c('headlinesink', class(hls))
  return(hls)
}

#'
#' @param linesink line-sink analytic element of class `linesink` or inherits from it.
#'
#' @return complex potential influence of `linesink` evaluated at points `x y`.
#' @noRd
#'
omegainf.linesink <- function(linesink, x, y, ...) {
  zeta <- x + y * 1i
  Z <- (2 * zeta - (linesink$z0 + linesink$z1)) / (linesink$z1 - linesink$z0)
  tol <- 1e-12
  zp1 <- ifelse(abs(Z + 1) < tol, tol, Z + 1)
  zm1 <- ifelse(abs(Z - 1) < tol, tol, Z - 1)
  omi <- linesink$L / (4*pi) * (zp1*log(zp1) - zm1*log(zm1))
  return(omi)
}

#'
#' @param linesink line-sink analytic element of class `linesink` or inherits from it.
#'
#' @return complex discharge influence of `linesink` evaluated at points `x y`.
#' @noRd
#'
domegainf.linesink <- function(linesink, x, y, ...) {
  zeta <- x + y * 1i
  Z <- (2 * zeta - (linesink$z0 + linesink$z1)) / (linesink$z1 - linesink$z0)
  tol <- 1e-12
  zp1 <- ifelse(abs(Z + 1) < tol, tol, Z + 1)
  zm1 <- ifelse(abs(Z - 1) < tol, tol, Z - 1)
  wi <- -(linesink$L / (2*pi)) * ((log(zp1) - log(zm1)) / (linesink$z1 - linesink$z0))
  return(wi)
}
