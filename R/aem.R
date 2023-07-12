
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
#' @return Object of class `aem` which is a list consisting of `k`, `top`, `base`, `n`,
#'    a list containing all elements, and a logical `solved` indicating if the model is solved.
#' @details If an element of class `headequation` is supplied, [solve.aem()] is called on the `aem`
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
#' @param x numeric x coordinates to evaluate `omega` at
#' @param y numeric y coordinates to evaluate `omega` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @export
#' @rdname omega
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
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `potential` at
#' @param y numeric y coordinates to evaluate `potential` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @export
#' @rdname potential
#' @include equation.R
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
#'
#' potential(ml, c(50, 0), c(25, -25))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' potential(ml, xg, yg, as.grid = TRUE)
#'
potential.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  pt <- Re(omega(aem, x, y, as.grid = as.grid, ...))
  return(pt)
}

#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `streamfunction` at
#' @param y numeric y coordinates to evaluate `streamfunction` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @export
#' @rdname streamfunction
#' @include equation.R
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
#'
#' streamfunction(ml, c(50, 0), c(25, -25))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' streamfunction(ml, xg, yg, as.grid = TRUE)
#'
streamfunction.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  sf <- Im(omega(aem, x, y, as.grid = as.grid, ...))
  return(sf)
}

#' Calculate the hydraulic head
#'
#' [heads()] computes the hydraulic head at the given x and y coordinates.
#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `heads` at
#' @param y numeric y coordinates to evaluate `heads` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @details Not to be confused with [utils::head()], which returns the first part of an object.
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the hydraulic head values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions c(`length(y), length(x)`) described by
#'     marginal vectors `x` and `y` containing the hydraulic head values at the grid points.
#' @export
#'
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
#'
#' heads(ml, c(50, 0), c(25, -25))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' heads(ml, xg, yg, as.grid = TRUE)
#'
heads <- function(aem, x, y, as.grid = FALSE, ...) { # rename as heads ???
  # TODO implement unconfined/confined flow
  hd <- potential(aem, x, y, as.grid = as.grid, ...) / (aem$k * (aem$top - aem$base))
  return(hd)
}

#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `domega` at
#' @param y numeric y coordinates to evaluate `domega` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @export
#' @rdname domega
#' @include equation.R
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
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


# TODO Qz
#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate `discharge` at
#' @param y numeric y coordinates to evaluate `discharge` at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param magnitude logical, should the magnitude of the discharge vector be returned as well? Default to `FALSE`. See details.
#' @param ... ignored
#'
#' @export
#' @rdname discharge
#' @include equation.R
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, uf)
#'
#' discharge(ml, c(50, 0), c(25, -25), magnitude = TRUE)
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#' discharge(ml, xg, yg, as.grid = TRUE)
#'
discharge.aem <- function(aem, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  W <- domega(aem, x, y, as.grid = as.grid, ...)
  Qx <- Re(W)
  Qy <- -Im(W)
  if(magnitude) {
    qv <- c(Qx, Qy, sqrt(Qx^2 + Qy^2))
    ndim <- 3
    nms <- c('Qx', 'Qy', 'Q')
  } else {
    qv <- c(Qx, Qy)
    ndim <- 2
    nms <- c('Qx', 'Qy')
  }
  if(as.grid) {
    Q <- array(qv, dim = c(dim(Qx), ndim), dimnames = list(NULL, NULL, nms)) # as used by {image} or {contour}. NROW and NCOL are switched
    Q <- image_to_matrix(Q)
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}

#' Title
#'
#' @param aem
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
darcy <- function(aem, x, y, ...) {
  # TODO add z
  Q <- discharge(aem, x, y, ...)
  b <- satthick(aem, x, y, ...)
  q <- Q / b
  return(q)
}

#' Title
#'
#' @param aem
#' @param x
#' @param y
#' @param n
#' @param b
#' @param R
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
velocity <- function(aem, x, y, R = 1, ...) {
  # TODO add z
  q <- darcy(aem, x, y, ...)
  v <- q / (aem$n * R)
  return(v)
}

#' Title
#'
#' @param aem
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
satthick <- function(aem, x, y, as.grid = FALSE, ...) {
  # TODO adjust for unconfined flow
  # TODO as.grid = TRUE
  d <- aem$top - aem$base
  b <- rep(d, max(length(x), length(y)))
  return(b)
}

#' Solve an `aem` model
#'
#' [solve.aem()] solves system of equations as constructed by the supplied  elements in the `aem` model
#'
#' @param aem `aem` object
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
solve.aem <- function(aem, ...) {
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
