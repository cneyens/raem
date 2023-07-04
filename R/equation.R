
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
  rhs <- element$pc
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

#' Calculate the complex potential
#'
#' [omega()] computes the complex potential for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex potential values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex potential values at the grid points.
#' @export
#' @rdname omega
#' @seealso [head.aem()], [potential()], [streamfunction()], [disc()], [disvec()]
#'
omega <- function(...) UseMethod('omega')

#' Calculate the discharge potential
#'
#' [potential()] computes the discharge potential for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the discharge potential values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the discharge potential values at the grid points.
#' @export
#' @rdname potential
#' @seealso [head.aem()], [omega()], [streamfunction()], [disc()], [disvec()]
#'
potential <- function(...) UseMethod('potential')

#' Calculate the streamfunction
#'
#' [streamfunction()] computes the streamfunction for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the streamfunction values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the streamfunction values at the grid points.
#' @export
#' @rdname streamfunction
#' @seealso [head.aem()], [omega()], [potential()], [disc()], [disvec()]
#'
streamfunction <- function(...) UseMethod('streamfunction')

#' Calculate the complex discharge
#'
#' [disc()] computes the complex discharge for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex discharge values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex discharge values at the grid points.
#' @export
#' @rdname disc
#' @seealso [head.aem()], [omega()], [potential()], [streamfunction()], [disvec()]
#'
disc <- function(...) UseMethod('disc')

#' Calculate the discharge vector
#'
#' [disvec()] computes the `Qx, Qy` and `Qz` components of the discharge vector for an `aem` or `element` object
#'     at the given x and y coordinates.
#'
#' @param ... ignored
#'
#' @return A matrix of dimensions `c(length(x), 3)` with the `Qx, Qy` and `Qz` components of
#'     the discharge vector at `x` and `y`. If `as.grid = TRUE`, an array of dimensions
#'     `c(length(y), length(x), 3)` described by marginal vectors `x` and `y` (columns and rows)
#'     containing the `Qx, Qy` and `Qz` components of the discharge vector.
#'
#' If `magnitude = TRUE`, the last dimension of the returned matrix is expanded to include
#'     the magnitude of the discharge vector, calculated as `sqrt(Qx^2 + Qy^2 + Qz^2)`.
#'
#' @export
#' @rdname disvec
#' @seealso [head.aem()], [omega()], [potential()], [streamfunction()], [disvec()]
#'
disvec <- function(...) UseMethod('disvec')

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
#' @seealso [omegainf()], [discinf()]
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
#' @seealso [potinf()], [discinf()]
#'
omegainf <- function(...) UseMethod('omegainf')

#' Calculate the complex discharge influence
#'
#' [discinf()] computes the complex discharge influence at the given x and y coordinates.
#' @param x numeric x coordinates to evaluate `discinf` at
#' @param y numeric y coordinates to evaluate `discinf` at
#' @param ... ignored
#'
#' @return A vector of `length(x)` (equal to `length(y)`) with the complex discharge influence values at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the complex discharge influence values at the grid points.
#' @noRd
#' @seealso [omegainf()], [potinf()]
#'
discinf <- function(...) UseMethod('omegainf')
