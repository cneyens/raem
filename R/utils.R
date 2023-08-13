
#' Compute the saturated thickness
#'
#' [satthick()] computes the saturated thickness of the aquifer of an `aem` object
#'     at the given x and y coordinates.
#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate at
#' @param y numeric y coordinates to evaluate at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @details If the aquifer is confined at `x` and `y`, the saturated thickness equals the aquifer thickness.
#'    For flow with variable saturated thickness, if the aquifer is unconfined at `x` and `y`, the saturated thickness
#'    it is calculated as the hydraulic head at `x` and `y` minus the aquifer base.
#'
#' @return [satthick()] returns a vector of `length(x)` (equal to `length(y)`) with the saturated thicknesses at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the saturated thicknesses at the grid points.
#' @export
#' @seealso [flow()], [state-variables()]
#' @examples
#' uf <- uniformflow(100, 0.001, 0)
#' rf <- constant(-1000, 0, 11)
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, uf, rf, type = 'confined')
#'
#' satthick(m, x = c(-200, 0, 200), y = 0) # confined
#' satthick(m, x = seq(-500, 500, length = 100),
#'          y = seq(-250, 250, length = 100), as.grid = TRUE)
#'
#' mv <- aem(k = 10, top = 10, base = 0, n = 0.2, uf, rf, type = 'variable')
#' satthick(mv, x = c(-200, 0, 200), y = 0) # variable
#'
satthick <- function(aem, x, y, as.grid = FALSE, ...) {

  if(as.grid) {
    df <- expand.grid(x = x, y = y) # increasing x and y values
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }

  if(aem$type == 'confined') {
    d <- aem$top - aem$base
  } else if(aem$type == 'variable') {
    h <- heads(aem, x = gx, y = gy)
    d <- ifelse(h >= aem$top, aem$top - aem$base, h - aem$base)
  }
  mb <- c(cbind(x = gx, y = gy, b = d)[,'b'], use.names = FALSE) # recycle x and y
  if(as.grid) {
    mb <- matrix(mb, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    mb <- image_to_matrix(mb)
  }
  return(mb)
}


#' Compute flow at a point in the direction of a given angle
#'
#' [dirflow()] computes a flow variable at the given points in the direction of the supplied angle
#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `flow` at
#' @param y numeric y coordinates to evaluate `flow` at
#' @param angle numeric, angle of the direction to evaluate `flow`, in degrees counterclockwise from the x-axis
#' @param flow character specifying which flow variable to use. Possible values are `discharge` (default), `darcy` and `velocity`. See [flow()].
#' @param as.grid logical, should a matrix be returned? Defaults to FALSE. See details.
#' @param ... additional arguments passed to [discharge()], [darcy()] or [velocity()]
#'
#' @details The x and y components of `flow` are used to calculate the directed value using `angle`.
#'    The `z` coordinate in [discharge()], [darcy()] or [velocity()] is set at the aquifer base. Under Dupuit-Forchheimer,
#'    the x and y components of flow do not change along the vertical axis.
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the flow values at `x` and `y` in the direction of `angle`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the directed flow values at the grid points.
#' @export
#' @seealso [flow()], [flow_through_line()]
#'
#' @examples
#' rf <- constant(-1000, 0, hc = 10)
#' uf <- uniformflow(TR = 100, gradient = 0.001, angle = -45)
#' w <- well(10, -50, Q = 200)
#'
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, uf)
#' dirflow(m, x = c(0, 100), y = 50, angle = -45)
#'
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, uf, w, type = 'confined')
#' dirflow(m, x = c(0, 50, 100), y = c(0, 50), angle = -90, flow = 'velocity', as.grid = TRUE)
#'
dirflow <- function(aem, x, y, angle,
                     flow = c('discharge', 'darcy', 'velocity'), as.grid = FALSE, ...) {
  theta <- angle * pi/180
  z <- aem$base
  flow <- match.arg(flow)
  flux <- switch(flow,
                 discharge = discharge(aem, x, y, z, as.grid = as.grid, ...),
                 darcy = darcy(aem, x, y, z, as.grid = as.grid, ...),
                 velocity = velocity(aem, x, y, z, as.grid = as.grid, ...))

  if(as.grid) {
    ndim <- length(dim(flux))
    if(ndim == 4) {
      f <- flux[,,,1] * cos(theta) + flux[,,,2] * sin(theta)
    } else if(ndim == 3) {
      f <- flux[,,1] * cos(theta) + flux[,,2] * sin(theta)
    }
  } else {
    f <- flux[,1] * cos(theta) + flux[,2] * sin(theta)
  }

  return(unname(f))
}

#' Calculate the total flow passing through a line
#'
#' [flow_through_line()] computes the integrated flow passing through a straight line at a right angle.
#'
#' @param aem `aem` object
#' @param x0 numeric, starting x location of line
#' @param y0 numeric, starting y location of line
#' @param x1 numeric, ending x location of line
#' @param y1 numeric, ending y location of line
#' @param flow character specifying which flow variable to use. Possible values are `discharge` (default) and `darcy`. See [flow()].
#' @param split logical, should the flow be split up into positive and negative flows or should they be summed (default)? See details.
#' @param ... ignored
#'
#' @details The flow is computed normal to the line and integrated along the line length using [stats::integrate()].
#'    The flow value is positive going to the left when looking in the direction of the line (i.e. to the left going from `x0-y0`
#'    to `x1-y1`).
#'
#' If `split = FALSE` (the default), a single value is returned which is the sum of the positive and negative flows perpendicular to the line.
#'    If `split = TRUE`, both the positive and negative component of the total flow through the line are returned.
#'
#' @return If `split = FALSE`, a single value with the total flow of variable `flow` passing through the line at a right angle.
#'    If `split = TRUE` a named vector with the total positive and total negative value of `flow` passing through the line.
#' @export
#' @seealso [flow()], [dirflow()]
#' @importFrom stats integrate
#' @examples
#' rf <- constant(-1000, 0, hc = 10)
#' uf <- uniformflow(TR = 100, gradient = 0.001, angle = -45)
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, uf)
#'
#' xg <- seq(-500, 500, l=100); yg <- seq(-300, 300, l=100)
#' contours(m, xg, yg, col='dodgerblue3', nlevels=20)
#'
#' x0 <- -200
#' y0 <- -50
#' x1 <- 300
#' y1 <- 100
#' lines(matrix(c(x0, y0, x1, y1), ncol = 2, byrow = TRUE))
#'
#' flow_through_line(m, x0, y0, x1, y1)
#' flow_through_line(m, x1, y1, x0, y0) # reverse direction of line
#'
#' w <- well(125, 200, 150)
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, uf, w)
#' contours(m, xg, yg, col='dodgerblue3', nlevels=20)
#' lines(matrix(c(x0, y0, x1, y1), ncol = 2, byrow = TRUE))
#'
#' flow_through_line(m, x0, y0, x1, y1, flow = 'darcy')
#' flow_through_line(m, x0, y0, x1, y1, flow = 'darcy', split = TRUE)
#'
flow_through_line <- function(aem, x0, y0, x1, y1, flow = c('discharge', 'darcy'), split = FALSE, ...) {

  # TODO allow matrix for multiple lines (vectorize integral and sum total flow)
  flow <- match.arg(flow) # no velocity allowed
  theta_line <- atan2(y1 - y0, x1 - x0)
  angle_norm <- (theta_line * 180/pi) + 90 # + 90: (= -90 * -flux) means positive flux is to the left when looking in direction of line
  len <- sqrt((y1 - y0)^2 + (x1 - x0)^2)
  nflow_int <- function(l, dropneg) {
    xi <- l * cos(theta_line) + x0
    yi <- l * sin(theta_line) + y0
    f <- dirflow(aem, xi, yi, angle = angle_norm, flow = flow)
    if(dropneg) f <- ifelse(f < 0, 0, f)
    return(f)
  }
  if(!split) {
    fint <- integrate(nflow_int, 0, len, dropneg = FALSE)$value
  } else {
    fintg <- integrate(nflow_int, 0, len, dropneg = FALSE)$value
    fintp <- integrate(nflow_int, 0, len, dropneg = TRUE)$value
    fintn <- fintg - fintp
    fint <- c('positive' = fintp, 'negative' = fintn)
  }
  return(fint)
}
