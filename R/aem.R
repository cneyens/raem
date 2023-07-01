
aem <- function(TR, ...) {
  l <- list(...)
  # if(length(l) == 1 && inherits(l[[1]], 'list')) l <- l[[1]]
  # https://stackoverflow.com/questions/72894519/deparse-and-substitute-on-ellipsis-to-get-names-of-parameters
  names(l) <- sapply(substitute(...()), deparse)
  if(any(vapply(l, function(i) !inherits(i, 'element'), FALSE))) stop('All supplied objects should be analytic elements', call. = FALSE)
  aem <- list(TR = TR, elements = l, solved = FALSE)
  class(aem) <- 'aem'

  if(any(vapply(aem$elements, function(i) inherits(i, 'headequation'), TRUE))) {
    aem <- solve(aem)
  }

  return(aem)
}

omega <- function(...) UseMethod('omega')
potential <- function(...) UseMethod('potential')
streamfunction <- function(...) UseMethod('streamfunction')
potinf <- function(...) UseMethod('potinf')
omegainf <- function(...) UseMethod('omegainf')

# complex discharge W at given point
# superposition of complex discharges of all elements
# complex discharge functions
disc <- function(...) UseMethod('disc')
# disc.element # not necessary


# discharge vector
# Qx = Re(W); Qy = -Im(W)
# returns matrix with two columns (Qx, Qy) when as.grid = FALSE
# otherwise 3D matrix with x,y,2 dimensions ([,,1] = Qx, [,,2] = Qy)
# no {discharge} functions per element TYPE (e.g. 'uf' or 'well'): they dispatch on element class
# magnitude = TRUE returns an additional column or layer with the magnitude of the vector computed as
# Q = sqrt(Qx^2 + Qy^2)
disvec <- function(...) UseMethod('disvec')

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
    om <- t(om[,dim(om)[2]:1])
  }
  return(om)
}
potential.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  pt <- Re(omega(aem, x, y, as.grid = as.grid))
  return(pt)
}
streamfunction.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  sf <- Im(omega(aem, x, y, as.grid = as.grid))
  return(sf)
}
head.aem <- function(aem, x, y, as.grid = FALSE, ...) { # rename as heads ???
  hd <- potential(aem, x, y, as.grid = as.grid) / aem$TR
  return(hd)
}

disc.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  w <- lapply(aem$elements, disc, x = gx, y = gy)
  w <- colSums(do.call(rbind, w))
  # w <- 0 + 0i
  # for(i in aem$elements) w <- w + disc(i, gx, gy)

  if(as.grid) {
    w <- matrix(w, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    w <- t(w[,dim(w)[2]:1])
  }
  return(w)
}


disvec.aem <- function(aem, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  W <- disc(aem, x, y, as.grid = as.grid, ...)
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
    Q <- aperm(Q[,dim(Q)[2]:1,], c(2,1,3))
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}



solve.aem <- function(aem, ...) {
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
