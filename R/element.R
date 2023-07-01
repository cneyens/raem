
element <- function(p, un = 0, ...) {
  el <- list()
  el$parameter <- p
  el$nunknowns <- un
  class(el) <- c('element')
  return(el)
}

omega.element <- function(element, x, y, ...) {
  om <- element$parameter * omegainf(element, x, y)
  return(om)
}

potential.element <- function(element, x, y, ...) {
  pt <- Re(omega(element, x, y))
  return(pt)
}

potinf.element <- function(element, x, y, ...) {
  pti <- Re(omegainf(element, x, y))
  return(pti)
}

disvec.element <- function(element, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  W <- disc(element, x, y, ...)
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
    Q <- aperm(Q[,dim(Q)[2]:1,], c(2,1,3))
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}

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
