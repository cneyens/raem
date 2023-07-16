
#' Create a analytic element model
#'
#' [aem()] creates an analytic element model to which elements can be added
#'
#' @param k numeric, hydraulic conductivity of aquifer
#' @param top numeric, top elevation of aquifer
#' @param base numeric, base elevation of aquifer
#' @param n numeric, effective porosity of aquifer as a fraction of total unit volume. Used for determining flow velocities with [velocity()].
#' @param ... objects of class `element`, or a single list with `element` objects
#'
#' @return [aem()] returns an object of class `aem` which is a list consisting of `k`, `top`, `base`, `n`,
#'    a list containing all elements, and a logical `solved` indicating if the model is solved.
#' @details When calling [aem()], if an element of class `headequation` is supplied, [solve.aem()] is called on the `aem`
#'     object before it is returned.
#' @export
#'
#' @examples
#' k <- 10
#' top <- 10
#' base <- 0
#' n <- 0.2
#' TR <- k * (top - base)
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' rf <- constant(xc = -500, yc = 0, h = 20)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
#'
#' aem(k, top, base, n, w, rf, uf)
#' aem(k, top, base, n, list(w, rf, uf))
#'
aem <- function(k, top, base, n, ...) {
  # TODO remove TR from list
  l <- list(...)
  if(length(l) == 1 && inherits(l[[1]], 'list') && !inherits(l[[1]], 'element')) l <- l[[1]]
  names(l) <- sapply(substitute(...()), deparse)
  if(any(vapply(l, function(i) !inherits(i, 'element'), FALSE))) stop('All supplied elements should be of class {element}', call. = FALSE)
  if(inherits(k, 'element') || !is.numeric(k)) stop('k should be numeric, not of class {element}')
  if(inherits(top, 'element') || !is.numeric(top)) stop('top should be numeric, not of class {element}')
  if(inherits(base, 'element') || !is.numeric(base)) stop('base should be numeric, not of class {element}')
  if(inherits(n, 'element') || !is.numeric(n)) stop('n should be numeric, not of class {element}')

  aem <- list(k = k, top = top, base = base, n = n, elements = l, solved = FALSE)
  class(aem) <- 'aem'

  if(any(vapply(aem$elements, function(i) inherits(i, 'headequation'), TRUE))) {
    aem <- solve(aem)
  }

  return(aem)
}

#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate at
#' @param y numeric y coordinates to evaluate at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @export
#' @rdname state-variables
#' @name state-variables
#' @include equation.R
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
#'
#' omega(ml, c(50, 0), c(25, -25))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' omega(ml, xg, yg, as.grid = TRUE)
#'
omega.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y) # increasing x and y values
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  om <- lapply(aem$elements, omega, x = gx, y = gy)
  om <- colSums(do.call(rbind, om))
  # om <- 0 + 0i
  # for(i in aem$elements) om <- om + omega(i, gx, gy)

  if(as.grid) {
    om <- matrix(om, nrow = length(x), ncol = length(y)) # as used by {image} or {contour}. NROW and NCOL are switched
    om <- image_to_matrix(om)
  }
  return(om)
}

#'
#' @export
#' @rdname state-variables
#' @name state-variables
#' @include equation.R
#' @examples
#' potential(ml, c(50, 0), c(25, -25))
#' potential(ml, xg, yg, as.grid = TRUE)
#'
potential.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  pt <- Re(omega(aem, x, y, as.grid = as.grid, ...))
  return(pt)
}

#'
#' @export
#' @rdname state-variables
#' @name state-variables
#' @include equation.R
#' @examples
#' streamfunction(ml, c(50, 0), c(25, -25))
#' streamfunction(ml, xg, yg, as.grid = TRUE)
#'
streamfunction.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  sf <- Im(omega(aem, x, y, as.grid = as.grid, ...))
  return(sf)
}

#'
#' @description [heads()] computes the hydraulic head at the given x and y coordinates.
#'
#' @details [heads()] should not to be confused with [utils::head()], which returns the first part of an object.
#' @return For [heads()], the same as for [omega()] but containing the hydraulic head values
#'    evaluated at `x` and `y`, which are computed from [potential()] and the aquifer parameters using [potential_to_head()].
#' @export
#' @rdname state-variables
#' @name state-variables
#' @examples
#' heads(ml, c(50, 0), c(25, -25))
#' heads(ml, xg, yg, as.grid = TRUE)
#' # do not confuse heads() with utils::head, which will give error
#' \dontrun{
#' head(ml, c(50, 0), c(25, -25))
#' }
#'
heads <- function(aem, x, y, as.grid = FALSE, ...) {
  # TODO implement unconfined/confined flow
  # TODO use [potential_to_head()]
  hd <- potential(aem, x, y, as.grid = as.grid, ...) / (aem$k * (aem$top - aem$base))
  return(hd)
}

#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate at
#' @param y numeric y coordinates to evaluate at
#' @param as.grid logical, should a matrix be returned? Defaults to `FALSE`. See details.
#' @param ... ignored or arguments passed from [velocity()] or [darcy()] to [discharge()].
#'
#' @export
#' @rdname flow
#' @name flow
#' @include equation.R
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
#' @name flow
#' @include equation.R
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
  # TODO leakage at bottom
  # TODO unconfined flow
  ntotal <- vapply(aem$elements, function(i) ifelse(inherits(i, 'areasink'), i$parameter, 0), 0.0)
  Qz <- -(gz - aem$base) * sum(ntotal) # TODO unconfined flow

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
#' @export
#' @rdname flow
#' @name flow
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
#' @export
#' @rdname flow
#' @name flow
#' @examples
#' velocity(ml, c(50, 0), c(25, -25), c(10, 5), magnitude = TRUE, R = 5)
#' velocity(ml, xg, yg, 5, as.grid = TRUE, R = 5)
#'
velocity.aem <- function(aem, x, y, z, as.grid = FALSE, magnitude = FALSE, R = 1, ...) {
  q <- darcy(aem, x, y, z, as.grid = as.grid, magnitude = magnitude, ...)
  v <- q / (aem$n * R)
  return(v)
}

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
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, uf, rf)
#'
#' satthick(m, x = c(-200, 0, 200), y = 0) # confined
#' satthick(m, x = seq(-500, 500, length = 100),
#'          y = seq(-250, 250, length = 100), as.grid = TRUE)
#' # TODO add unconfined example
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

  # TODO adjust for unconfined flow
  d <- aem$top - aem$base
  mb <- cbind(x = gx, y = gy, b = d)[,'b'] # recycle x and y
  if(as.grid) {
    mb <- matrix(mb, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    mb <- image_to_matrix(mb)
  }
  return(mb)
}

#' Solve an `aem` model
#'
#' [solve.aem()] solves system of equations as constructed by the supplied  elements in the `aem` model
#'
#' @param a `aem` object
#' @param b ignored
#' @param ... ignored
#'
#' @details In order to be solved, `aem` should contain at least
#'     1 `headequation` element.
#'     [solve.aem()] sets up the system of equations, and calls [solve()] to
#'     solve.
#'     Constructing an `aem` object by a call to [aem()] automatically calls [solve.aem()].
#'
#' @return The solved `aem` object, i.e. after finding the solution
#'     to the system of equations as constructed by the contained elements.
#' @export
#'
#' @examples
#'
#' rf <- constant(-500, 0, 20)
#' hdw <- headwell(xw = 0, yw = 100, rw = 0.3, hw = 8)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2) |>
#'            add_element(rf, name = 'rf') |>
#'            add_element(hdw, name = 'headwell')
#' ml <- solve(ml)
#'
solve.aem <- function(a, b, ...) {
  aem <- a
  # TODO if no headequations, just return model instead of error
  if(!any(vapply(aem$elements, function(i) inherits(i, 'headequation'), TRUE))) {
    stop('Model should contain at least 1 headequation element in order to be solved', call. = FALSE)
  }

  nun <- vapply(aem$elements, function(i) i$nunknowns, 1)
  esolve <- aem$elements[which(nun == 1)]
  nunknowns <- length(esolve)
  m <- matrix(0, nrow = nunknowns, ncol = nunknowns)
  rhs <- rep(0, nunknowns)
  for(irow in 1:nunknowns) {
    eq <- equation(esolve[[irow]], aem)
    m[irow,] <- eq[[1]]
    rhs[irow] <- eq[[2]]
  }
  solution <- solve(m, rhs)
  for(irow in 1:nunknowns) esolve[[irow]]$parameter <- solution[irow]
  aem$elements[which(nun == 1)] <- esolve
  aem$solved <- TRUE
  return(aem)
}
