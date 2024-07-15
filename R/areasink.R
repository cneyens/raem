
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
#' @seealso [headareasink()]
#'
#' @examples
#' areasink(xc = -500, yc = 0, N = 0.001, R = 500)
#'
#' # flux assuming a constant head difference over a confining unit
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

#' Create a head-specified area-sink analytic element
#'
#' [headareasink()] creates a area-sink analytic element with constant specified head. The constant leakage flux
#'    into or out of the aquifer from the area-sink is computed by solving the corresponding `aem` model.
#'
#' @param xc numeric, x location of the center of the area-sink
#' @param yc numeric, y location of the center of the area-sink
#' @param hc numeric, specified hydraulic head at the center of the area-sink
#' @param R numeric, radius of the circular area-sink
#' @param resistance numeric, hydraulic resistance of the area-sink at its connection with the aquifer. Defaults to 0 (no resistance).
#' @param location character, either `top` (default) or `base` specifying the vertical position of the area-sink.
#' @param ... ignored
#'
#' @details The constant leakage flux from the area-sink is computed by solving the `aem` model given
#'    the specified head `hc` for the area-sink. This head is located at the so-called collocation point,
#'    which is placed at the center of the area-sink. A positive flux is into the aquifer. Note that this head-dependent
#'    flux is constant over the domain and computed only at the collocation point. The flux is therefore determined
#'    by the difference in aquifer head and specified head at that location only, and does not vary across the domain with
#'    varying aquifer head.
#'
#' The resistance can be increased for a area-sink in poor connection with the aquifer, e.g. because of
#'    a confining unit of low hydraulic conductivity between the aquifer and the area-sink. If the aquifer is unconfined
#'    (i.e. has a variable saturated thickness), the system of equations becomes non-linear with respect to the hydraulic head
#'    and iteration is required to solve the model.
#'
#' @return Circular head-specified area-sink analytic element which is an object of class `headareasink` and inherits from `areasink`.
#' @export
#' @seealso [areasink()]
#' @examples
#'
#' # using headareasink, the head difference depends on the computed instead of given aquifer head
#' headareasink(xc = -500, yc = 0, hc = 3, R = 500, res = 1000)
#' headareasink(xc = -500, yc = 0, hc = 3, R = 500, location = 'base')
#'
headareasink <- function(xc,
                         yc,
                         hc,
                         R,
                         resistance = 0,
                         location = c('top', 'base'),
                         ...) {
  has <- areasink(xc = xc, yc = yc, N = 0, R = R, location = location)
  has$hc <- hc
  has$nunknowns <- 1
  has$resistance <- resistance
  class(has) <- c('headareasink', class(has))
  return(has)
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
  Qr <- ifelse(r <= Rs, r/2, Rs^2/(2*r))
  W <- Qr*(x - areasink$xc)/r - 1i*(Qr*(y - areasink$yc)/r) # polar to cartesian
  W[r < 1e-15] <- 0 + 0i # continuous across boundary
  return(W)
}
