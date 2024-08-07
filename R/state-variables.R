
#' @title Calculate state-variables
#'
#' @description [heads()] computes the hydraulic head at the given x and y coordinates for an `aem` object.
#'
#' @param na.below logical indicating if calculated head values below the aquifer base should be set to `NA`. Defaults to `TRUE`. See [potential_to_head()].
#' @details [heads()] should not to be confused with [utils::head()], which returns the first part of an object.
#' @return For [heads()], a vector of `length(x)` (equal to `length(y)`) with the hydraulic head values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the hydraulic head values at the grid points.
#'     The heads are computed from [potential()] and the aquifer parameters using [potential_to_head()].
#' @export
#' @rdname state-variables
#' @name state-variables
#' @seealso [flow()], [satthick()], [head_to_potential()]
#' @examples
#' w <- well(xw = 55, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' rf <- constant(xc = -1000, yc = 1000, hc = 10)
#' ml <- aem(k = 10, top = 10, base = -15, n = 0.2, w, uf, rf)
#'
#' xg <- seq(-100, 100, length = 5)
#' yg <- seq(-75, 75, length = 3)
#'
#' # Hydraulic heads
#' heads(ml, c(50, 0), c(25, -25))
#' heads(ml, xg, yg, as.grid = TRUE)
#'
#' # do not confuse heads() with utils::head, which will give an error
#' try(
#' head(ml, c(50, 0), c(25, -25))
#' )
#'
heads <- function(aem, x, y, as.grid = FALSE, na.below = TRUE, ...) {
  phi <- potential(aem, x, y, as.grid = as.grid, ...)
  hd <- potential_to_head(aem, phi, na.below = na.below)
  return(hd)
}

#'
#' @description [omega()] computes the complex potential for an `aem` or `element` object
#'   at the given x and y coordinates.
#'
#' @return For [omega()], the same as for [heads()] but containing the complex potential values
#'   evaluated at `x` and `y`.
#' @export
#' @rdname state-variables
#'
omega <- function(...) UseMethod('omega')

#'
#' @description [potential()] computes the discharge potential for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [potential()], the same as for [heads()] but containing the discharge potential values
#'    evaluated at `x` and `y`, which are the real components of [omega()].
#' @export
#' @rdname state-variables
#'
potential <- function(...) UseMethod('potential')

#'
#' @description [streamfunction()] computes the stream function for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [streamfunction()], the same as for [heads()] but containing the stream function values
#'    evaluated at `x` and `y`, which are the imaginary components of [omega()].
#' @export
#' @rdname state-variables
#'
streamfunction <- function(...) UseMethod('streamfunction')

#' Calculate the potential influence
#'
#' [potinf()] computes the potential influence at the given x and y coordinates.
#'
#' @param ... ignored
#' @param x numeric x coordinates to evaluate `potinf` at.
#' @param y numeric y coordinates to evaluate `potinf` at.
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the potential influence values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the potential influence values at the grid points.
#' @noRd
#' @seealso [omegainf()], [domegainf()]
#'
potinf <- function(...) UseMethod('potinf')

#' Calculate the complex potential influence
#'
#' [omegainf()] computes the complex potential influence at the given x and y coordinates.
#'
#' @param ... ignored
#' @param x numeric x coordinates to evaluate `omegainf` at.
#' @param y numeric y coordinates to evaluate `omegainf` at.
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex potential influence values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex potential influence values at the grid points.
#' @noRd
#' @seealso [potinf()], [domegainf()]
#'
omegainf <- function(...) UseMethod('omegainf')

#'
#' @param aem `aem` object.
#' @param x numeric x coordinates to evaluate the variable at.
#' @param y numeric y coordinates to evaluate the variable at.
#' @param as.grid logical, should a matrix be returned? Defaults to `FALSE`. See details.
#' @param ... ignored
#'
#' @export
#' @rdname state-variables
#' @examples
#' # Complex potential
#' omega(ml, c(50, 0), c(25, -25))
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
#' @examples
#' # Discharge potential
#' potential(ml, c(50, 0), c(25, -25))
#'
potential.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  pt <- Re(omega(aem, x, y, as.grid = as.grid, ...))
  return(pt)
}

#'
#' @export
#' @rdname state-variables
#' @examples
#' # Stream function
#' streamfunction(ml, c(50, 0), c(25, -25))
#'
streamfunction.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  sf <- Im(omega(aem, x, y, as.grid = as.grid, ...))
  return(sf)
}


#'
#' @param element analytic element of class `element`.
#'
#' @export
#' @rdname state-variables
#' @examples
#' # For elements
#' omega(w, c(50, 0), c(-25, 25))
#'
omega.element <- function(element, x, y, ...) {
  if(length(element$parameter) == 1) {
    om <- element$parameter * omegainf(element, x, y, ...)
  } else { # sum of sn*omegainf_n, only for elements with more than 1 parameter
    om <- omegainf(element, x, y, ...) %*% element$parameter
  }

  return(c(om))
}

#'
#' @export
#' @rdname state-variables
#' @examples
#' potential(w, c(50, 0), c(-25, 25))
#'
potential.element <- function(element, x, y, ...) {
  pt <- Re(omega(element, x, y, ...))
  return(pt)
}

#'
#' @export
#' @rdname state-variables
#' @examples
#' streamfunction(w, c(50, 0), c(-25, 25))
#'
streamfunction.element <- function(element, x, y, ...) {
  sf <- Im(omega(element, x, y, ...))
  return(sf)
}

#'
#' @param element analytic element of class `element`.
#' @noRd
#'
potinf.element <- function(element, x, y, ...) {
  pti <- Re(omegainf(element, x, y, ...))
  return(pti)
}

#' Convert hydraulic head to potential and vice versa
#'
#' [head_to_potential()] calculates the discharge potential from the hydraulic head.
#'
#' @param aem `aem` object.
#' @param h numeric hydraulic head values as vector or matrix.
#' @param ... ignored
#'
#' @return [head_to_potential()] returns the discharge potentials calculated from `h`, in the same
#'    structure as `h`.
#' @export
#' @rdname head_to_potential
#' @examples
#' k <- 10
#' top <- 10; base <- 0
#' uf <- uniformflow(TR = 100, gradient = 0.001, angle = -45)
#' rf <- constant(TR, xc = -1000, yc = 0, hc = 10)
#' w1 <- well(200, 50, Q = 250)
#' m <- aem(k, top, base, n = 0.2, uf, rf, w1, type = 'variable') # variable saturated thickness
#' mc <- aem(k, top, base, n = 0.2, uf, rf, w1, type = 'confined') # constant saturated thickness
#' xg <- seq(-500, 500, length = 100)
#' yg <- seq(-250, 250, length = 100)
#'
#' h <- heads(m, x = xg, y = yg, as.grid = TRUE)
#' hc <- heads(mc, x = xg, y = yg, as.grid = TRUE)
#' pot <- head_to_potential(m, h)
#' potc <- head_to_potential(mc, hc)
#'
head_to_potential <- function(aem, h, ...) {
  b <- aem$top - aem$base
  if(aem$type == 'confined') {
    phi <- h * aem$k * b
  } else if(aem$type == 'variable') {
    cn <- 0.5 * aem$k * b^2 + aem$k * b * aem$base
    phi <- ifelse(h >= aem$top, h * aem$k * b - cn, 0.5 * aem$k * (h - aem$base)^2)
  }
  return(phi)
}

#' @description [potential_to_head()] calculates the hydraulic head from the discharge potential.
#'
#' @param phi numeric discharge potential values as vector or matrix.
#' @param na.below logical indicating if calculated head values below the aquifer base should be set to `NA`. Defaults to `TRUE`. See details.
#' @export
#' @return [potential_to_head()] returns the hydraulic heads calculated from `phi`, in the same
#'    structure as `phi`.
#'
#' The conversion of potential to head or vice versa is different for confined (constant saturated thickness)
#'    and unconfined (variable saturated thickness) aquifers as set by the `type` argument in `aem()`.
#'
#' If `na.below = FALSE`, negative potentials can be converted to hydraulic heads if flow is unconfined (`aem$type = 'variable'`).
#'    The resulting heads are below the aquifer base. This may be useful for some use cases, e.g. in preliminary model construction
#'    or for internal functions. In most cases however, these values should be set to `NA` (the default behavior) since other analytic
#'    elements will continue to extract or inject water even though the saturated thickness of the aquifer is negative,
#'    which is not realistic. In those cases, setting `aem$type = 'confined'` might prove useful. Also note that these heads below the
#'    aquifer base will not be correctly re-converted to potentials using [head_to_potential()]. As such, caution should be taken when
#'    setting `na.below = FALSE`.
#'
#' @rdname head_to_potential
#' @examples
#' phi <- potential(m, x = xg, y = yg, as.grid = TRUE)
#' phic <- potential(mc, x = xg, y = yg, as.grid = TRUE)
#' hds <- potential_to_head(m, phi)
#' hdsc <- potential_to_head(mc, phic)
#'
#' # Converting negative potentials results in NA's with warning
#' try(
#' potential_to_head(m, -300)
#' )
#'
#' # unless na.below = FALSE
#' potential_to_head(m, -300, na.below = FALSE)
#'
potential_to_head <- function(aem, phi, na.below = TRUE, ...) {
  b <- aem$top - aem$base
  if(aem$type == 'confined') {
    h <- phi / (aem$k * b)
  } else if(aem$type == 'variable') {
    phit <- 0.5 * aem$k * b^2
    cn <- phit + aem$k * b * aem$base
    if(na.below) {
      h <- ifelse(phi >= phit, (phi + cn)/(aem$k * b), sqrt(2*phi/aem$k) + aem$base)
    } else {
      h <- ifelse(phi >= phit, (phi + cn)/(aem$k * b), sqrt(2*abs(phi)/aem$k) * ifelse(phi < 0, 1i, 1))
      h <- ifelse(phi >= phit, h, ifelse(phi < 0, Im(Conj(h)), Re(h))  + aem$base)
    }
  }
  return(h)
}
