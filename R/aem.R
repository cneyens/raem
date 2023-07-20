
#' Create an analytic element model
#'
#' [aem()] creates an analytic element model to which elements can be added
#'
#' @param k numeric, hydraulic conductivity of aquifer
#' @param top numeric, top elevation of aquifer
#' @param base numeric, base elevation of aquifer
#' @param n numeric, effective porosity of aquifer as a fraction of total unit volume. Used for determining flow velocities with [velocity()].
#' @param ... objects of class `element`, or a single (named) list with `element` objects
#' @param type character specifying the type of flow in the aquifer, either `variable` (default) or `confined`. See details.
#'
#' @return [aem()] returns an object of class `aem` which is a list consisting of `k`, `top`, `base`, `n`,
#'    a list containing all elements with the names of the objects specified in `...`, and a logical `solved`
#'    indicating if the model is solved.
#' @details The default `type = 'variable'` allows for unconfined/confined flow, i.e. flow with variable saturated thickness. If `type = 'confined'`,
#'    the saturated thickness is always constant and equal to the aquifer thickness.
#'
#' When calling [aem()], if an element of class `headequation` is supplied, [solve.aem()] is called on the `aem`
#'     object before it is returned.
#' @export
#' @seealso [add_element()] [contours()]
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
#' hdw <- headwell(xw = 0, yw = 100, rw = 0.3, hc = 8)
#' ls <- linesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, sigma = 1)
#'
#' # Creating aem ----
#' m <- aem(k, top, base, n, w, rf, uf, hdw, ls)
#'
#' # or with elements in named list
#' m <- aem(k, top, base, n,
#'          list('well' = w, 'constant' = rf, 'flow' = uf, 'headwell' = hdw, 'river' = ls),
#'          type = 'confined')
#'
aem <- function(k, top, base, n, ..., type = c('variable', 'confined')) {

  type <- match.arg(type)

  l <- list(...)
  # setting names
  nms <- sapply(substitute(...()), deparse) # object names
  lnames <- names(l) # named arguments
  if(length(lnames) == 0) lnames <- nms
  noname <- which(nchar(lnames) == 0) # substitue lnames with nms
  if(length(noname) > 0) lnames[noname] <- nms[noname]
  names(l) <- lnames

  if(length(l) == 1 && inherits(l[[1]], 'list') && !inherits(l[[1]], 'element')) {
    l <- l[[1]]
  }

  # checks
  if(any(vapply(l, function(i) !inherits(i, 'element'), FALSE))) stop('All supplied elements should be of class \'element\'', call. = FALSE)
  if(inherits(k, 'element') || !is.numeric(k)) stop('k should be numeric, not of class \'element\'', call. = FALSE)
  if(inherits(top, 'element') || !is.numeric(top)) stop('top should be numeric, not of class \'element\'', call. = FALSE)
  if(inherits(base, 'element') || !is.numeric(base)) stop('base should be numeric, not of class \'element\'', call. = FALSE)
  if(inherits(n, 'element') || !is.numeric(n)) stop('n should be numeric, not of class \'element\'', call. = FALSE)
  # if(length(l) == 0) warning('No elements supplied. Add them using \'add_element\'', call. = FALSE)
  if(any(duplicated(names(l)))) stop('Duplicate names in \'elements\' not allowed', call. = FALSE)

  aem <- list(k = k, top = top, base = base, n = n, elements = l, type = type, solved = FALSE)
  class(aem) <- 'aem'

  aem <- solve(aem)

  return(aem)
}

#' Solve an `aem` model
#'
#' [solve.aem()] solves system of equations as constructed by the supplied  elements in the `aem` model
#'
#' @param a `aem` object
#' @param b ignored
#' @param ... ignored
#'
#' @details [solve.aem()] sets up the system of equations, and calls [solve()] to
#'     solve. If head-specified elements are supplied, an element of class `constant` as
#'     created by [constant()] (also called the reference point), should be supplied as well.
#'
#' Constructing an `aem` object by a call to [aem()] automatically calls [solve.aem()].
#'
#' @return The solved `aem` object, i.e. after finding the solution
#'     to the system of equations as constructed by the contained elements.
#' @export
#' @rdname aem
#' @examples
#' # Solving ----
#' m <- solve(m)
#'
#' # solving requires a reference point (constant) element if head-specified elements are supplied
#' try(
#'   m <- aem(k = k, top = top, base = base, n = n, w, uf, hdw)
#' )
#'
solve.aem <- function(a, b, ...) {
  aem <- a

  # no unknowns
  if(!any(vapply(aem$elements, function(i) i$nunknowns > 0, TRUE))) {
    aem$solved <- TRUE
    return(aem)
  }
  # check if reference point is provided
  # this is actually only necessary when headlinesink elements are specified
  if(!any(vapply(aem$element, function(i) inherits(i, 'constant'), TRUE))) {
    stop('Please provide an element of class \'constant\' when solving with head-specified elements',
         call. = FALSE)
  }

  nun <- vapply(aem$elements, function(i) i$nunknowns, 1)
  esolve <- aem$elements[which(nun > 0)]
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
  aem$elements[which(nun > 0)] <- esolve
  aem$solved <- TRUE
  return(aem)
}

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
#'     to be solved by [solve.aem()].
#'
#' @return a list containing the value of the stiffness matrix (A) and the RHS (b)
#'     for a given element.
#' @noRd
#' @seealso [solve.aem()]
#'
equation <- function(element, aem, ...) {
  if(!(element$nunknowns > 0)) stop('element should have 1 or more unknowns', call. = FALSE)
  row <- vector(mode = 'numeric')
  rhs <- head_to_potential(aem, element$hc)
  xc <- element$xc
  yc <- element$yc
  for(i in aem$elements) {
    if(i$nunknowns > 0) {
      row[length(row)+1] <- potinf(i, xc, yc)
    } else {
      rhs <- rhs - potential(i, xc, yc)
    }
  }
  return(list(row, rhs))

}

#' Create a base element
#'
#' [element()] creates a base element with a parameter and number of unknowns.
#'
#' @param p numeric parameter value
#' @param un numeric number of unknowns
#' @param ... ignored
#'
#' @return An object of class `element` which is a list with `p` and `un` elements.
#'     This is used in constructing analytic elements by adding other variabels to this base class.
#' @noRd
#'
element <- function(p, un = 0, ...) {
  el <- list()
  el$parameter <- p
  el$nunknowns <- un
  class(el) <- c('element')
  return(el)
}

#' Add element to existing `aem` object
#'
#' @param aem `aem` object
#' @param element analytic element of class `element`
#' @param name optional name of the element as character. Duplicate names in `aem` are not allowed.
#' @param solve logical, should the model be solved after adding the new element? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @return The `aem` model with the addition of `element`. If `solve = TRUE`, the model is solved using [solve.aem()].
#' @export
#' @seealso [aem()]
#' @examples
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2)
#' add_element(m, constant(xc = 0, yc = 1000, hc = 12), name = 'rf')
#' @examplesIf R.version$major >= 4 & R.version$minor >= 1
#' # add_element() is pipe-friendly
#' aem(k = 10, top = 10, base = 0, n = 0.2) |>
#'     add_element(constant(xc = 0, yc = 1000, hc = 12),
#'                 name = 'rf') |>
#'     add_element(headwell(xw = 0, yw = 100, rw = 0.3, hc = 8),
#'                 name = 'headwell', solve = TRUE)
#'
#'
add_element <- function(aem, element, name = NULL, solve = FALSE, ...) {
  if(!inherits(aem, 'aem')) stop('\'aem\' object should be of class aem', call. = FALSE)
  if(!inherits(element, 'element')) stop('\'element\' object should be of class element', call. = FALSE)
  if(is.null(aem$elements)) aem$elements <- list()
  aem$elements[[length(aem$elements) + 1]] <- element
  if(!is.null(name)) {
    if(name %in% names(aem$elements)) stop('element ', '\'', name, '\'', ' already exists', call. = FALSE)
    names(aem$elements)[length(aem$elements)] <- name
  }
  if(solve) aem <- solve(aem)
  return(aem)
}
