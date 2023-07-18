
#' @title Calculate flow variables
#'
#' @description [domega()] computes the complex discharge for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [domega()],  vector of `length(x)` (equal to `length(y)`) with the complex discharge values at `x` and `y`,
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex discharge values at the grid points.
#'     [domega()] is the derivative of [omega()] in the x and y directions.
#' @rdname flow
#' @name flow
#' @export
#' @seealso [state-variables()], [satthick()], [head_to_potential()]
#'
domega <- function(...) UseMethod('domega')

#'
#' @description [discharge()] computes the `x, y` and `z` components of the discharge vector for an `aem` or `element` object
#'     at the given x, y and z coordinates.
#'
#' @details There is no [discharge()], [darcy()] or [velocity()] method for an object of class `element` because an `aem` object is required
#'    to obtain the aquifer base and top.
#'
#' @return For [discharge()], a matrix with the number of rows equal to the number of points to evaluate
#'    the discharge vector at, and with columns `Qx`, `Qy` and `Qz` corresponding to `x, y` and `z` components
#'    of the discharge vector at coordinates `x`, `y` and `z`. If `as.grid = TRUE`, an array of dimensions
#'    `c(length(y), length(x), length(z), 3)` described by marginal vectors `x`, `y` and `z` (columns, rows and third dimension)
#'    containing the `x, y` and `z` components of the discharge vector (`Qx`, `Qy` and `Qz`) as the fourth dimension.
#'
#' The `x` component of [discharge()] is the real value of [domega()], the `y` component
#'    the imaginary component and the `z` component is calculated based on area-sink strengths.
#'
#' If `magnitude = TRUE`, the last dimension of the returned array is expanded to include
#'     the magnitude of the discharge/darcy/velocity vector, calculated as `sqrt(Qx^2 + Qy^2 + Qz^2)`
#'     (or `sqrt(qx^2 + qy^2 + qz^2)` or `sqrt(vx^2 + vy^2 + vz^2)`, respectively).
#' @rdname flow
#' @export
#'
discharge <- function(...) UseMethod('discharge')

#'
#' @description [darcy()] computes the `x, y` and `z` components of the Darcy flux vector (also called specific discharge vector)
#'    for an `aem` object at the given x, y and z coordinates.
#'
#' @export
#' @return For [darcy()], the same as for [discharge()] but with the `x`, `y` and `z` components of the
#'    Darcy flux vector (`qx`, `qy` and `qz`). The value are computed by dividing the values of [discharge()] by
#'    the saturated thickness at `x`, `y` and `z`.
#' @rdname flow
#'
darcy <- function(...) UseMethod('darcy')

#'
#' @description [velocity()] computes the `x, y` and `z` components of the average linear groundwater flow velocity vector
#'     for an `aem` object at the given x, y and z coordinates.
#' @rdname flow
#' @export
#' @return For [velocity()], the same as for [discharge()] but with the `x`, `y` and `z` components of the
#'    average linear groundwater flow velocity vector (`vx`, `vy` and `vz`). The values are computed by dividing
#'    the [darcy()] values by the effective porosity (`aem$n`) and the retardation coefficient `R`.
#'
velocity <- function(...) UseMethod('velocity')

#' Calculate the complex discharge influence
#'
#' [domegainf()] computes the complex discharge influence at the given x and y coordinates.
#' @param x numeric x coordinates to evaluate `domegainf` at
#' @param y numeric y coordinates to evaluate `domegainf` at
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex discharge influence values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex discharge influence values at the grid points.
#' @noRd
#' @seealso [omegainf()], [potinf()]
#'
domegainf <- function(...) UseMethod('domegainf')

#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate at
#' @param y numeric y coordinates to evaluate at
#' @param as.grid logical, should a matrix be returned? Defaults to `FALSE`. See details.
#' @param ... ignored or arguments passed from [velocity()] or [darcy()] to [discharge()].
#'
#' @rdname flow
#' @export
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' as <- areasink(xc = 0, yc = 0, N = 0.001, R = 500)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf, as)
#'
#' domega(ml, c(50, 0), c(25, -25))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' domega(ml, xg, yg, as.grid = TRUE)
#'
domega.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  w <- lapply(aem$elements, domega, x = gx, y = gy)
  w <- colSums(do.call(rbind, w))
  # w <- 0 + 0i
  # for(i in aem$elements) w <- w + domega(i, gx, gy)

  if(as.grid) {
    w <- matrix(w, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    w <- image_to_matrix(w)
  }
  return(w)
}

#'
#' @param z numeric z coordinates to evaluate at
#' @param magnitude logical, should the magnitude of the flow vector be returned as well? Default to `FALSE`. See details.
#' @param verbose logical, if `TRUE` (default), warnings with regards to setting `Qz` to NA are printed. See details.
#'
#' @details If the `z` coordinate is above the saturated aquifer level (i.e. the water-table for unconfined conditions or
#'    the aquifer top for confined conditions), or below the aquifer base, `Qz` values are set to NA with a warning (if `verbose = TRUE`).
#'    The `Qx` and `Qy` values are not set to NA, for convenience in specifying the `z` coordinate when only lateral flow
#'    is of interest.
#' @export
#' @rdname flow
#' @examples
#' discharge(ml, c(50, 0), c(25, -25), z = ml$top)
#' discharge(ml, c(50, 0), c(25, -25), z = c(ml$top, 5), magnitude = TRUE)
#' discharge(ml, xg, yg, z = ml$top, as.grid = TRUE)
#' discharge(ml, c(50, 0), c(25, -25), z = ml$top + c(0, 0.5)) # NA for z > top
#'
discharge.aem <- function(aem, x, y, z, as.grid = FALSE, magnitude = FALSE, verbose = TRUE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y, z = z)
    gx <- df$x
    gy <- df$y
    gz <- df$z
  } else {
    gx <- x
    gy <- y
    gz <- z
  }

  # Get Qx and Qy
  W <- domega(aem, gx, gy, as.grid = FALSE, ...)
  Qx <- Re(W)
  Qy <- -Im(W)

  # Get Qz: depends on area-sinks
  # TODO unconfined flow
  ntop <- vapply(aem$elements, function(i) ifelse(inherits(i, 'areasink') && i$location == 'top', i$parameter, 0), 0.0)
  nbase <- vapply(aem$elements, function(i) ifelse(inherits(i, 'areasink') && i$location == 'base', i$parameter, 0), 0.0)
  sat <- satthick(aem, gx, gy)

  Qz <- (gz - aem$base) * (-sum(ntop) - sum(nbase)) + sum(nbase)*sat

  # set Qz to NA if z coordinate above saturated part or below aquifer base
  # TODO keep this behaviour? Shouldn't Qx and Qy also be set to NA?
  outside_v <- outside_vertical(aem, gx, gy, gz)$outside
  if(any(outside_v) && verbose) {
    warning('Setting Qz values to NA for z above saturated aquifer level or below aquifer base', call. = FALSE)
  }
  Qz <- ifelse(outside_v, NA, Qz)
  # Qx <- ifelse(outside_v, NA, Qx)
  # Qy <- ifelse(outside_v, NA, Qy)

  if(magnitude) {
    qv <- cbind(Qx, Qy, Qz, sqrt(Qx^2 + Qy^2 + Qz^2))
    ndim <- 4
    nms <- c('Qx', 'Qy', 'Qz', 'Q')
  } else {
    qv <- cbind(Qx, Qy, Qz)
    ndim <- 3
    nms <- c('Qx', 'Qy', 'Qz')
  }
  if(as.grid) {
    Q <- array(c(qv), dim = c(length(x), length(y), length(z), ndim), dimnames = list(NULL, NULL, NULL, nms)) # as used by {image} or {contour}. NROW and NCOL are switched
    Q <- image_to_matrix(Q)
  } else {
    Q <- matrix(c(qv), ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}

#'
#' @rdname flow
#' @export
#' @examples
#' darcy(ml, c(50, 0), c(25, -25), c(10, 5), magnitude = TRUE)
#' darcy(ml, xg, yg, 10, as.grid = TRUE)
#'
darcy.aem <- function(aem, x, y, z, as.grid = FALSE, magnitude = FALSE, ...) {
  Q <- discharge(aem, x, y, z, as.grid = as.grid, magnitude = magnitude, ...)
  b <- satthick(aem, x, y, as.grid = as.grid, ...)
  q <- Q / array(b, dim = dim(Q))
  return(q)
}

#' @param R numeric, retardation coefficient. Defaults to 1 (no retardation).
#'
#' @rdname flow
#' @export
#' @examples
#' velocity(ml, c(50, 0), c(25, -25), c(10, 5), magnitude = TRUE, R = 5)
#' velocity(ml, xg, yg, 5, as.grid = TRUE, R = 5)
#'
velocity.aem <- function(aem, x, y, z, as.grid = FALSE, magnitude = FALSE, R = 1, ...) {
  q <- darcy(aem, x, y, z, as.grid = as.grid, magnitude = magnitude, ...)
  v <- q / (aem$n * R)
  return(v)
}

#'
#' @param element analytic element of class `element`
#'
#' @export
#' @rdname flow
#' @examples
#' # For elements
#' domega(w, c(50, 0), c(-25, 25))
#'
domega.element <- function(element, x, y, ...) {
  wi <- element$parameter * domegainf(element, x, y, ...)
  return(wi)
}
