
#' Obtain parameters for linear system of equations
#'
#' [equation()] obtains the values of the stiffness matrix and the RHS vector
#'     for a given element used to construct the linear system of equations.
#'
#' @param element analytic element for which the coefficients should be determined.
#' @param aem `aem` object
#' @param ... ignored
#'
#' @details [equation()] is used to set up the linear system of equations `Ax = b`
#'     to be solved by [solve()].
#'
#' @return a list containing the value of the stiffness matrix (A) and the RHS (b)
#'     for a given element.
#' @noRd
#' @seealso [solve()]
#'
equation <- function(element, aem, ...) {
  # TODO return 0 if !inherits(element, 'headequation')
  if(!inherits(element, 'headequation')) stop('element should be of class headequation', call. = FALSE)
  row <- vector(mode = 'numeric')
  rhs <- head_to_potential(aem, element$hc)
  xc <- element$xc
  yc <- element$yc
  for(i in aem$elements) {
    if(i$nunknowns == 1) {
      row[length(row)+1] <- potinf(i, xc, yc)
    } else {
      rhs <- rhs - potential(i, xc, yc)
    }
  }
  return(list(row, rhs))

}

#' @title Calculate state-variables
#'
#' @description [omega()] computes the complex potential for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [omega()], a vector of `length(x)` (equal to `length(y)`) with the complex potential values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex potential values at the grid points.
#' @export
#' @rdname state-variables
#' @name state-variables
#' @seealso [flow()], [satthick()], [head_to_potential()]
#'
omega <- function(...) UseMethod('omega')

#'
#' @description [potential()] computes the discharge potential for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [potential()], the same as for [omega()] but containing the discharge potential values
#'    evaluated at `x` and `y`, which are the real components of [omega()].
#' @export
#' @rdname state-variables
#' @name state-variables
#'
potential <- function(...) UseMethod('potential')

#'
#' @description [streamfunction()] computes the streamfunction for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [streamfunction()], the same as for [omega()] but containing the streamfunction values
#'    evaluated at `x` and `y`, which are the imaginary components of [omega()].
#' @export
#' @rdname state-variables
#' @name state-variables
#'
streamfunction <- function(...) UseMethod('streamfunction')

#' @title Calculate flow variables
#'
#' @description [domega()] computes the complex discharge for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @return For [domega()],  vector of `length(x)` (equal to `length(y)`) with the complex discharge values at `x` and `y`,
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex discharge values at the grid points.
#'     [domega()] is the derivate of [omega()] in the x and y directions.
#'
#' @export
#' @rdname flow
#' @name flow
#' @seealso [state-variables()], [satthick()], [head_to_potential()]
#'
domega <- function(...) UseMethod('domega')

#'
#' @description [discharge()] computes the `x, y` and `z` components of the discharge vector for an `aem` or `element` object
#'     at the given x, y and z coordinates.
#'
#' @return For [discharge()], a matrix of dimensions `c(length(x), 3)` with the `x, y` and `z` components (`Qx`, `Qy` and `Qz`)
#'    of the discharge vector at coordinates `x`, `y` and `z`. If `as.grid = TRUE`, an array of dimensions
#'    `c(length(y), length(x), 3)` described by marginal vectors `x` and `y` (columns and rows)
#'    containing the `x, y` and `z` components of the discharge vector (`Qx`, `Qy` and `Qz`).
#'
#' The `x` component of [discharge()] is the real value of [domega()], the `y` component
#'    the imaginary component and `z` is calculated based on area-sink strengths.
#'
#' If `magnitude = TRUE`, the last dimension of the returned matrix is expanded to include
#'     the magnitude of the discharge/darcy/velocity vector, calculated as `sqrt(Qx^2 + Qy^2 + Qz^2)`
#'     (or `sqrt(qx^2 + qy^2 + qz^2)` or `sqrt(vx^2 + vy^2 + vz^2)`, respectively).
#'
#' @export
#' @rdname flow
#' @name flow
#'
discharge <- function(...) UseMethod('discharge')

#'
#' @description [darcy()] computes the `x, y` and `z` components of the Darcy flux vector (also called specific discharge vector)
#'    for an `aem` object at the given x, y and z coordinates.
#'
#' @details There is no [darcy()] or [velocity()] method for an object of class `element` because an `aem` object is required
#'    to calculate the saturated thickness using [satthick()].
#' @export
#' @return For [darcy()], the same as for [discharge()] but with the `x`, `y` and `z` components of the
#'    Darcy flux vector (`qx`, `qy` and `qz`). The value are computed by dividing the values of [discharge()] by the saturated thickness.
#' @rdname flow
#' @name flow
#'
darcy <- function(...) UseMethod('darcy')

#'
#' @description [velocity()] computes the `x, y` and `z` components of the average linear groundwater flow velocity vector
#'     for an `aem` object at the given x, y and z coordinates.
#'
#' @export
#' @return For [velocity()], the same as for [discharge()] but with the `x`, `y` and `z` components of the
#'    average linear groundwater flow velocity vector (`vx`, `vy` and `vz`). The values are computed by dividing
#'    the [darcy()] values by the effective porosity and the retardation coefficient.
#' @rdname flow
#' @name flow
#'
velocity <- function(...) UseMethod('velocity')

#' Calculate the potential influence
#'
#' [potinf()] computes the potential influence at the given x and y coordinates.
#'
#' @param ... ignored
#' @param x numeric x coordinates to evaluate `potinf` at
#' @param y numeric y coordinates to evaluate `potinf` at
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
#' @param x numeric x coordinates to evaluate `omegainf` at
#' @param y numeric y coordinates to evaluate `omegainf` at
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex potential influence values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex potential influence values at the grid points.
#' @noRd
#' @seealso [potinf()], [domegainf()]
#'
omegainf <- function(...) UseMethod('omegainf')

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

#' Title
#'
#' @param aem
#' @param h
#' @param ...
#'
#' @return
#' @export
#' @rdname head_to_potential
#' @examples
head_to_potential <- function(aem, h, ...) {
  # TODO unconfined flow
  # export ??
  pot <- h * aem$k * (aem$top - aem$base)
  return(pot)
}

#' Title
#'
#' @param aem
#' @param pot
#' @param ...
#'
#' @return
#' @export
#' @rdname head_to_potential
#' @examples
potential_to_head <- function(aem, pot, ...) {
  # TODO unconfined flow
  h <- pot / (aem$k * (aem$top - aem$base))
  return(h)
}
