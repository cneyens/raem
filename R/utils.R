
#' Compute the saturated thickness
#'
#' [satthick()] computes the saturated thickness of the aquifer from an `aem` object
#'     at the given x and y coordinates.
#'
#' @param aem `aem` object.
#' @param x numeric x coordinates to evaluate at.
#' @param y numeric y coordinates to evaluate at.
#' @param as.grid logical, should a matrix be returned? Defaults to `FALSE`. See details.
#' @param ... additional arguments passed to [heads()] when `aem$type = 'variable'`.
#'
#' @details If the aquifer is confined at `x` and `y`, the saturated thickness equals the aquifer thickness.
#'    For flow with variable saturated thickness (`aem$type = 'variable'`), if the aquifer is unconfined at `x` and `y`,
#'    the saturated thickness is calculated as the hydraulic head at `x` and `y` minus the aquifer base.
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the saturated thicknesses at `x` and `y`.
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
#' s <- satthick(m, x = seq(-500, 500, length = 100),
#'               y = seq(-250, 250, length = 100), as.grid = TRUE)
#' str(s)
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
    h <- heads(aem, x = gx, y = gy, as.grid = FALSE, ...)
    d <- ifelse(h >= aem$top, aem$top - aem$base, h - aem$base)
  }
  mb <- c(cbind(x = gx, y = gy, b = d)[,'b'], use.names = FALSE) # recycle x and y
  if(as.grid) {
    mb <- matrix(mb, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    mb <- image_to_matrix(mb)
  }
  return(mb)
}


#' Compute flow in the direction of a given angle
#'
#' [dirflow()] computes a flow variable at the given points in the direction of the supplied angle.
#'
#' @param aem `aem` object.
#' @param x numeric x coordinates to evaluate `flow` at.
#' @param y numeric y coordinates to evaluate `flow` at.
#' @param angle numeric, angle of the direction to evaluate `flow`, in degrees counterclockwise from the x-axis.
#' @param flow character specifying which flow variable to use. Possible values are `discharge` (default), `darcy` and `velocity`. See [flow()].
#' @param as.grid logical, should a matrix be returned? Defaults to FALSE. See details.
#' @param ... additional arguments passed to [discharge()], [darcy()] or [velocity()].
#'
#' @details The x and y components of `flow` are used to calculate the directed value using `angle`.
#'    The `z` coordinate in [discharge()], [darcy()] or [velocity()] is set at the aquifer base. Under Dupuit-Forchheimer,
#'    the x and y components of the flow vector do not change along the vertical axis.
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
#' dirflow(m, x = c(0, 50, 100), y = c(0, 50), angle = -90,
#' flow = 'velocity', as.grid = TRUE)
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
#' @param x0 numeric, starting x location of line.
#' @param y0 numeric, starting y location of line.
#' @param x1 numeric, ending x location of line.
#' @param y1 numeric, ending y location of line.
#' @param flow character specifying which flow variable to use. Possible values are `discharge` (default) and `darcy`. See [flow()].
#' @param split logical, should the flow be split up into positive and negative flows (`TRUE`) or should they be summed (`FALSE`; default)? See details.
#' @param ... ignored
#'
#' @details The flow is computed normal to the line and integrated along the line length using [stats::integrate()].
#'    The flow value is positive going to the left when looking in the direction of the line (i.e. to the left going from `x0-y0`
#'    to `x1-y1`).
#'
#' If `split = FALSE` (the default), a single value is returned which is the sum of the positive and negative flows perpendicular to the line.
#'    If `split = TRUE`, both the positive and negative component of the total flow through the line are returned.
#'
#' If the line corresponds to a line element, the integration might fail. Try to perturbate the line vertices slightly in that case.
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
#' contours(m, xg, yg, col='dodgerblue', nlevels=20)
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
#' contours(m, xg, yg, col='dodgerblue', nlevels=20)
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

#' Get the computed discharge from an element
#'
#' [element_discharge()] obtains the computed discharge into or out of the aquifer
#'    for a individual analytic element or all elements of a given type.
#'
#' @param aem `aem` object.
#' @param name character vector with the name of the element(s) as available in `aem$elements`.
#' @param type character with the type (class) of element to obtain the summed discharge from. See details.
#' @param ... ignored
#'
#' @details Either `name` or `type` should be specified. If `type` is specified, only one type is allowed.
#'    Possible values are `'headwell', 'well', 'linesink', 'headlinesink', 'areasink'` or `'headareasink'`.
#'
#' Only elements that add or remove water from the aquifer will return a non-zero discharge value.
#'
#' @return A numeric named vector of length `length(name)` with the discharge into (negative) or out of (positive)
#'    the aquifer. If `type` is specified, a single named numeric value with the total discharge into (negative) or
#'    out of (positive) the aquifer which is the sum of all individual elements of class `type`.
#' @export
#'
#' @examples
#' k <- 10
#' top <- 10
#' base <- 0
#' n <- 0.2
#' TR <- k * (top - base)
#'
#' rf <- constant(xc = -500, yc = 0, h = 20)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
#' w1 <- well(xw = 50, yw = 0, Q = 200)
#' w2 <- well(xw = 0, yw = 100, Q = 400)
#' hw <- headwell(xw = -100, yw = 0, hc = 7.5)
#' hls <- headlinesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, hc = 8)
#' as <- areasink(xc = 0, yc = 0, N = 0.0005, R = 500)
#' m <- aem(k, top, base, n, rf, uf, w1, w2, hw, hls, as)
#'
#' element_discharge(m, name = c('hls', 'as'))
#' element_discharge(m, type = 'well')
#'
#' # zero discharge for uniform flow element as it does not add or remove water
#' element_discharge(m, name = 'uf')
element_discharge <- function(aem, name = NULL, type = NULL, ...) {
  if((is.null(name) && is.null(type)) || (!is.null(name) && !is.null(type))) stop('Either "name" or "type" should be specified', call. = FALSE)
  if(!is.null(name)) {
    elnames <- names(aem$elements)
    id <- which(elnames %in% name)
    if(any(length(id) == 0)) stop('element with name "', name[which(length(id) == 0)], '" not found', call. = FALSE)
    get_Q_name <- function(el) {
      if(inherits(el, 'well')) {
        return(el$parameter)
      } else if(inherits(el, 'linesink')) {
        return(el$parameter * el$L)
      } else if(inherits(el, 'areasink')) {
        return(-el$parameter * pi * el$R^2)
      } else {
        return(0) # uniformflow, constant, linedoublet
      }
    }
    Q <- unlist(lapply(aem$elements[id], get_Q_name))
    names(Q) <- name
  } else {
    types <- c('headwell', 'well', 'linesink', 'headlinesink', 'areasink', 'headareasink')
    if(length(type) > 1) stop('Only one type should be specified', call. = FALSE)
    if(!(type %in% types)) stop('"type" should be one of ', paste(types, collapse = ', '), call. = FALSE)

    get_Q_type <- function(el) {
      if(inherits(el, type, which = TRUE) == 1) { # makes sure el inherits from type but not if it is a parent class
        if(grepl('well', type)) {
          return(el$parameter)
        } else if(grepl('linesink', type)) {
          return(el$parameter * el$L)
        } else if(grepl('areasink', type)) {
          return(-el$parameter * pi * el$R^2)
        }
      } else {
        return(0)
      }
    }
    Q <- sum(unlist(lapply(aem$elements, get_Q_type)))
    names(Q) <- type
  }
  return(Q)
}


#' Determine if point is inside polygon
#'
#' [point_in_polygon()] determines if a point is inside a polygon. Points on the line are treated
#'    as outside. This is an implementation of a ray-casting algorithm.
#'
#' @param point numeric vector of length 2 with x and y coordinates of point
#' @param polygon matrix with columns x and y containing the vertices of the polygon. It should not be closed.
#'
#' @return A logical indicating if the point lies inside the polygon.
#' @noRd
#'
#' @examples
#' poly <- rbind(c(-0.25, 0.25), c(0.32, 0.64), c(-0.25, 1.5), c(1, 1), c(0.76, 0.3),
#'               c(0.25, -0.25), c(-0.25, -0.5))
#'
#' pinside <- c(0.5,0.5) # point inside the poly
#' poutside <- c(1.5,1) # point outside the poly
#'
#' point_in_polygon(pinside, poly)
#' point_in_polygon(poutside, poly)
#'
#' n <- 1000
#' rpts <- cbind(x = runif(n, -0.5, 1.5), y = runif(n, -1, 2))
#' io <- apply(rpts, 1, point_in_polygon, polygon = poly)
#' df <- cbind(as.data.frame(rpts), col = ifelse(io, 'blue', 'red'))
#'
#' plot(df$x, df$y, col = df$col, xlim = c(-0.6, 1.6), ylim = c(-1.1, 2.1))
#' polygon(poly)
#'
point_in_polygon <- function(point, polygon) {
  x <- point[1]
  y <- point[2]
  poly_x <- polygon[,1]
  poly_y <- polygon[,2]
  nvert <- length(poly_x)
  inside <- FALSE
  j <- nvert
  for (i in 1:nvert) {
    if (((poly_y[i] > y) != (poly_y[j] > y)) &&
        (x < (poly_x[j] - poly_x[i]) * (y - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])) {
      inside <- !inside
    }
    j <- i
  }
  return(inside)
}
