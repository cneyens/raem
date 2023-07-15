
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

#'
#' @param element analytic element of class `element`
#'
#' @export
#' @rdname state-variables
#' @name state-variables
#' @examples
#' # For elements
#' omega(w, c(50, 0), c(-25, 25))
#'
omega.element <- function(element, x, y, ...) {
  om <- element$parameter * omegainf(element, x, y, ...)
  return(om)
}

#'
#' @export
#' @rdname state-variables
#' @name state-variables
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
#' @name state-variables
#' @examples
#' streamfunction(w, c(50, 0), c(-25, 25))
#'
streamfunction.element <- function(element, x, y, ...) {
  sf <- Im(omega(element, x, y, ...))
  return(sf)
}

#'
#' @param element analytic element of class `element`
#' @noRd
#'
potinf.element <- function(element, x, y, ...) {
  pti <- Re(omegainf(element, x, y, ...))
  return(pti)
}

#'
#' @param element analytic element of class `element`
#'
#' @export
#' @rdname flow
#' @name flow
#' @examples
#' # For elements
#' domega(w, c(50, 0), c(-25, 25))
#'
domega.element <- function(element, x, y, ...) {
  wi <- element$parameter * domegainf(element, x, y, ...)
  return(wi)
}

#'
#' @export
#' @rdname flow
#' @name flow
#' @examples
#' discharge(w, c(50, 0), c(-25, 25))
#' discharge(w, c(50, 0), c(-25, 25), as.grid = TRUE, magnitude = TRUE)
#'
discharge.element <- function(element, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  # TODO add z
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  W <- domega(element, x, y, ...)
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
    Q <- array(qv, dim = c(length(x), length(y), ndim), dimnames = list(NULL, NULL, nms)) # as used by {image} or {contour}. NROW and NCOL are switched
    Q <- image_to_matrix(Q)
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
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
#'
#' @examples
#' aem(k = 10, top = 10, base = 0, n = 0.2) |>
#'     add_element(constant(xc = 0, yc = 1000, hc = 12),
#'                 name = 'rf') |>
#'     add_element(headwell(xw = 0, yw = 100, rw = 0.3, hw = 8),
#'                 name = 'headwell', solve = TRUE)
#'
add_element <- function(aem, element, name = NULL, solve = FALSE, ...) {
  if(!inherits(aem, 'aem')) stop('{aem} should be of class aem', call. = FALSE)
  if(!inherits(element, 'element')) stop('{element} should be of class element', call. = FALSE)
  if(is.null(aem$elements)) aem$elements <- list()
  aem$elements[[length(aem$elements) + 1]] <- element
  if(!is.null(name)) {
    if(name %in% names(aem$elements)) stop('element ', '\'', name, '\'', ' already exists', call. = FALSE)
    names(aem$elements)[length(aem$elements)] <- name
  }
  if(inherits(element, 'headequation') & solve) aem <- solve(aem)
  return(aem)
}
