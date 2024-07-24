
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
#' @param verbose logical indicating if information during the solving process should be printed. Defaults to `FALSE`.
#' @param maxiter integer specifying the maximum allowed iterations for a non-linear solution. Defaults to 10. See [solve.aem()].
#'
#' @return [aem()] returns an object of class `aem` which is a list consisting of `k`, `top`, `base`, `n`,
#'    a list containing all elements with the names of the objects specified in `...`, and a logical `solved`
#'    indicating if the model is solved.
#' @details The default `type = 'variable'` allows for unconfined/confined flow, i.e. flow with variable saturated thickness. If `type = 'confined'`,
#'    the saturated thickness is always constant and equal to the aquifer thickness.
#'
#' [solve.aem()] is called on the `aem` object before it is returned, which solves the system of equations.
#'
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
aem <- function(k, top, base, n, ..., type = c('variable', 'confined'), verbose = FALSE, maxiter = 10) {

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

  aem <- solve(aem, maxiter = maxiter, verbose = verbose)

  return(aem)
}

#' Solve an `aem` model
#'
#' [solve.aem()] solves the system of equations as constructed by the elements in the `aem` model
#'
#' @param a `aem` object
#' @param b ignored
#' @param maxiter integer specifying the maximum allowed iterations for a non-linear solution. Defaults to 10. See details.
#' @param verbose logical indicating if information during the solving process should be printed. Defaults to `FALSE`.
#' @param ... ignored
#'
#' @details [solve.aem()] sets up the system of equations, and calls [solve()] to
#'     solve. If head-specified elements are supplied, an element of class `constant` as
#'     created by [constant()] (also called the reference point), should be supplied as well.
#'
#' Constructing an `aem` object by a call to [aem()] automatically calls [solve.aem()].
#'
#' If the system of equations is non-linear, for example when the flow system is unconfined (variable
#'    saturated thickness) and elements with hydraulic resistance are specified, a Picard iteration is entered.
#'    During each Picard iteration step (outer iteration), the previously solved model parameters are used to set up and
#'    solve a linear system of equations. The model parameters are then updated and the next outer iteration step is
#'    entered, until `maxiter` iterations are reached. For an linear model, `maxiter` is ignored.
#'
#' @return [solve.aem()] returns the solved `aem` object, i.e. after finding the solution
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
solve.aem <- function(a, b, maxiter = 10, verbose = FALSE, ...) {
  aem <- a
  if(length(aem$elements) == 0) {
    if(verbose) warning('No elements in model. Can not solve.', call. = FALSE)
    return(aem)
  }
  if(verbose) cat('Solving analytic element model ...', '\n')

  # no unknowns
  if(!any(vapply(aem$elements, function(i) i$nunknowns > 0, TRUE))) {
    if(verbose) cat(' Linear model with', length(aem$elements), ifelse(length(aem$elements) == 1, 'element', 'elements'),
                    'and 0 unknowns', '\n', 'Model solved', '\n')
    aem$solved <- TRUE
    return(aem)
  }

  # check if reference point is provided
  # this is actually only necessary when headlinesink elements are specified
  if(!any(vapply(aem$element, function(i) inherits(i, 'constant'), TRUE))) {
    stop('Please provide an element of class \'constant\' when solving for unknown parameters in elements',
         call. = FALSE)
  }

  # TODO allow for multiple unknowns per element (have to change the loops and omega())
  nun <- vapply(aem$elements, function(i) i$nunknowns, 1)
  esolve_id <- which(nun > 0)
  esolve <- aem$elements[esolve_id]
  nunknowns <- sum(nun)
  is_nonlinear <- any(vapply(esolve, function(i) ifelse(is.null(i$resistance), 0, i$resistance), 0) != 0) && aem$type == 'variable'
  if(!is_nonlinear) maxiter <- 1
  if(verbose) {
    cat(ifelse(is_nonlinear, ' Non-linear', ' Linear'), 'model with', length(aem$elements), ifelse(length(aem$elements) == 1, 'element', 'elements'),
        'and', nunknowns, ifelse(nunknowns == 1, 'unknown', 'unknowns'), '\n')
  }

  # TODO closer criterion to exit Picard loop when criterion is satisfied
  # e.g. if max absolute head difference at control points at iter i and iter i-1 < hclose, exit Picard loop

  # Picard iteration
  if(verbose & is_nonlinear) cat(' Entering outer iteration loop ...', '\n')
  for(iter in seq_len(maxiter)) {
    if(verbose & is_nonlinear) cat('  Iteration', iter, '\n')

    # set up system of equations
    m <- matrix(0, nrow = nunknowns, ncol = nunknowns)
    rhs <- rep(0, nunknowns)
    irow <- 0
    for(i in seq_along(esolve)) {
      el <- esolve[[i]]
      nunel <- el$nunknowns
      irow <- seq(1, nunel) + irow[length(irow)]
      eq <- equation(el, aem, esolve_id[i])
      m[irow,] <- eq[[1]]
      rhs[irow] <- eq[[2]]
    }

    # solve and set model parameters
    solution <- solve(m, rhs)
    irow <- 0
    for(i in seq_along(esolve)) {
      nunel <- esolve[[i]]$nunknowns
      irow <- seq(1, nunel) + irow[length(irow)]
      esolve[[i]]$parameter <- solution[irow]
    }
    aem$elements[esolve_id] <- esolve

  }

  if(verbose) cat('Model solved', '\n')
  aem$solved <- TRUE
  # aem$linear <- !is_nonlinear
  return(aem)
}

#' Obtain parameters for linear system of equations
#'
#' [equation()] obtains the values of the stiffness matrix and the RHS vector
#'     for a given element used to construct the linear system of equations.
#'
#' @param element analytic element for which the coefficients should be determined.
#' @param aem `aem` object
#' @param id integer indicating which place `element` is in `aem$elements`
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
equation <- function(element, aem, id, ...) {
  if(!(element$nunknowns > 0)) stop('element should have 1 or more unknowns', call. = FALSE)
  row <- vector(mode = 'numeric')
  xc <- element$xc
  yc <- element$yc
  resf <- resfac(element, aem)

  if(inherits(element, 'linedoublet')) {
    rhs <- 0
    tol <- 1e-12
    theta_norm <- atan2(Im(element$z1 - element$z0), Re(element$z1 - element$z0)) - pi/2
    xci <- xc - tol * cos(theta_norm)
    yci <- yc - tol * sin(theta_norm)
    xco <- xc + tol * cos(theta_norm)
    yco <- yc + tol * sin(theta_norm)

    for(i in aem$elements) {
      nun <- i$nunknowns
      if(nun > 0) {
        Qinf <- domegainf(i, xc, yc)
        dphi <- potinf(i, xci, yci) - potinf(i, xco, yco)
        row <- c(row, Re(Qinf)*cos(theta_norm) - Im(Qinf)*sin(theta_norm) - resf*dphi)

      } else {
        Q <- domega(i, xc, yc)
        dphi <- potential(i, xci, yci) - potential(i, xco, yco)
        rhs <- rhs - Re(Q)*cos(theta_norm) - Im(Q)*sin(theta_norm) + resf*dphi
      }
    }

  } else {
    if(inherits(element, 'inhomogeneity')) {
      rhs <- 0
    } else {
      rhs <- head_to_potential(aem, element$hc)
    }

    for(e in seq_along(aem$elements)) {
      i <- aem$element[[e]]
      nun <- i$nunknowns
      if(nun > 0) {
        if(e == id) {
          resfelm = resf
        } else {
          resfelm = 0
        }
        row <- c(row, potinf(i, xc, yc) - resfelm)
      } else {
        rhs <- rhs - potential(i, xc, yc)
      }
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

#' Get the resistance factor of an analytic element
#'
#' @param element analytic element with unknowns
#' @param aem `aem` object
#'
#' @return Numeric vector with the resistance factor at every collocation point of the element
#' @noRd
#'
resfac <- function(element, aem) {

  if(element$nunknowns == 0) stop('nunknowns should be > 0 to get resfac', call. = FALSE)

  if(inherits(element, 'inhomogeneity')) {
    resfac <- element$k / (element$k - aem$k) # TODO check this

  } else if(inherits(element, 'linedoublet')) {
    if(element$resistance == Inf) { # Impermeable wall
      resfac <- rep(0, length(element$xc))
    } else {
      if(aem$type == 'confined') {
        b <- satthick(aem, element$xc, element$yc)
        bsat <- b
      } else {
        tol <- 1e-12
        theta_norm <- atan2(Im(element$z1 - element$z0), Re(element$z1 - element$z0)) - pi/2
        xci <- element$xc - tol * cos(theta_norm)
        yci <- element$yc - tol * sin(theta_norm)
        xco <- element$xc + tol * cos(theta_norm)
        yco <- element$yc + tol * sin(theta_norm)

        hi <- heads(aem, c(xci), c(yci), na.below = FALSE)
        ho <- heads(aem, c(xco), c(yco), na.below = FALSE)
        hi <- ifelse(hi >= aem$top, aem$top, hi)
        ho <- ifelse(ho >= aem$top, aem$top, ho)

        b <- 0.5*(hi + ho) - aem$base
        bsat <- satthick(aem, element$xc, element$yc, na.below = FALSE)
      }
      if(element$resistance == 0) element$resistance <- 1e-12
      resfac <- bsat / (element$resistance * b * aem$k)
    }


  } else if(inherits(element, 'headwell')) {
    if(aem$type == 'confined') {
      b <- satthick(aem, element$xc, element$yc)
      bsat <- b
    } else {
      h <- heads(aem, element$xc, element$yc, na.below = FALSE)
      b <- ifelse(h >= aem$top, aem$top, 0.5*(h + element$hc)) - aem$base
      bsat <- satthick(aem, element$xc, element$yc, na.below = FALSE)
    }
    resfac <- element$resistance * aem$k * b / (2 * pi * element$rw * bsat)

  } else if(inherits(element, 'headlinesink')) {
    # remove width if width = 0
    width <- ifelse(element$width == 0, 1, element$width)

    # Haitjema, 1995, eq. 5.64 & 5.66
    if(aem$type == 'confined') {
      b <- satthick(aem, element$xc, element$yc)
    } else {
      h <- heads(aem, element$xc, element$yc, na.below = FALSE)
      b <- ifelse(h >= aem$top, aem$top, 0.5*(h + element$hc)) - aem$base
    }
    resfac <- element$resistance * aem$k * b / width

  } else if(inherits(element, 'headareasink')) {
    # Haitjema, 1995, eq. 5.64 & 5.66
    if(aem$type == 'confined') {
      b <- satthick(aem, element$xc, element$yc)
    } else {
      h <- heads(aem, element$xc, element$yc, na.below = FALSE)
      b <- ifelse(h >= aem$top, aem$top, 0.5*(h + element$hc)) - aem$base
    }
    resfac <- -element$resistance * aem$k * b

  } else {
    resfac <- rep(0, length(element$xc))
  }

  resfac <- ifelse(is.na(resfac), 0, resfac) # if bsat or b in denominator = 0 in first iteration
  return(resfac)
}

#' Add or remove an element from an existing `aem` object
#'
#' [add_element()] adds a new element to or from an `aem` object.
#'
#' @param aem `aem` object
#' @param element analytic element of class `element`
#' @param name optional name of the element as character. Duplicate names in `aem` are not allowed.
#' @param solve logical, should the model be solved after adding or removing the element? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @return The `aem` model with the addition of `element` or with the removal of an element. If `solve = TRUE`,
#'    the model is solved using [solve.aem()]. The name of the new element is taken from the `name` argument,
#'    the object name or set to `element_1` with `1` being the index of the new element in the element list. See examples.
#' @export
#' @seealso [aem()]
#' @examples
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2)
#' mnew <- add_element(m, constant(xc = 0, yc = 1000, hc = 12), name = 'rf')
#'
#' # if name not supplied, tries to obtain it from object name
#' rf <- constant(xc = 0, yc = 1000, hc = 12)
#' mnew <- add_element(m, rf)
#'
#' # or else sets it sequentially from number of elements
#' mnew <- add_element(m, constant(xc = 0, yc = 1000, hc = 12))
#' @examplesIf getRversion() >= '4.1.0'
#' # add_element() is pipe-friendly
#' mnew <- aem(k = 10, top = 10, base = 0, n = 0.2) |>
#'     add_element(rf, name = 'rf') |>
#'     add_element(headwell(xw = 0, yw = 100, rw = 0.3, hc = 8),
#'                 name = 'headwell', solve = TRUE)
#'
add_element <- function(aem, element, name = NULL, solve = FALSE, ...) {
  if(!inherits(aem, 'aem')) stop('\'aem\' object should be of class aem', call. = FALSE)
  if(!inherits(element, 'element')) stop('\'element\' object should be of class element', call. = FALSE)
  if(is.null(aem$elements)) aem$elements <- list()
  aem$elements[[length(aem$elements) + 1]] <- element
  if(is.null(name)) {
    name <- sapply(substitute(element), deparse) # object name
    if(length(name) > 1) name <- NULL # weak test ...
  }
  if(!is.null(name)) {
    if(name %in% names(aem$elements)) stop('Element ', '\'', name, '\'', ' already exists', call. = FALSE)
    names(aem$elements)[length(aem$elements)] <- name
  } else {
    names(aem$elements)[length(aem$element)] <- paste('element', length(aem$elements), sep = '_')
  }
  if(solve) {
    aem <- solve(aem)
  } else if(aem$solved) {
    aem$solved <- FALSE
  }
  return(aem)
}

#'
#' @description [remove_element()] removes an element from the `aem` object based on its name or type.
#'
#' @param type class of the element(s) to remove. Either `name` or `type` should be specified in [remove_element()].
#'
#' @export
#' @rdname add_element
#'
#' @examples
#' # removing elements
#' mnew <- remove_element(mnew, name = 'rf')
#' mnew <- remove_element(mnew, type = 'headwell')
#'
remove_element <- function(aem, name = NULL, type = NULL, solve = FALSE, ...) {
  if(!inherits(aem, 'aem')) stop('\'aem\' object should be of class aem', call. = FALSE)
  if((is.null(name) && is.null(type)) || (!is.null(name) && !is.null(type))) stop('Either "name" or "type" should be specified', call. = FALSE)

  if(is.null(aem$elements)) return(aem)

  if(is.null(name)) {
    istype <- vapply(aem$elements, function(i) inherits(i, type, which = TRUE), 0)
    id <- which(istype == 1)
    if(length(id) == 0) {
      warning('No elements of type "', type, '" found', call. = FALSE)
      return(aem)
    }
    aem$elements <- aem$elements[-id]
  } else {
    elnames <- names(aem$elements)
    id <- which(elnames %in% name)
    if(any(length(id) == 0)) {
      warning('Element with name "', name[which(length(id) == 0)], '" not found', call. = FALSE)
      return(aem)
    }
    aem$elements <- aem$elements[-id]
  }

  if(solve) {
    aem <- solve(aem)
  } else if(aem$solved) {
    aem$solved <- FALSE
  }
  return(aem)
}
