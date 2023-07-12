
#' Title
#'
#' @param well
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @noRd
#'
#' @examples
reached_well <- function(well, x, y, ...) {
  xw <- Re(well$zetaw)
  yw <- Im(well$zetaw)
  d <- sqrt((x - xw)^2 + (y - yw)^2) - well$rw
  # if(d <= 0) d <- 0
  return(d <= 0)
}

#' Title
#'
#' @param aem
#' @param x0
#' @param y0
#' @param times
#' @param forward
#' @param R
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
tracelines <- function(aem, x0, y0, times, forward = TRUE, R = 1, ...) {

  # wrapper to obtain velocity
  vxvy <- function(t, y, parms, ...) {
    m <- matrix(y, ncol = 2, byrow = TRUE) # necessary ??
    v <- velocity(parms$aem, x=m[,1], y=m[,2], R = parms$R)
    direction <- ifelse(forward, 1, -1)
    return(list(c(direction * v)))
  }

  # TODO other functions when particles should terminate
  rootfun <- function(t, y, parms) {
    m <- matrix(y, ncol = 2, byrow = TRUE) # necessary ??
    wls <- vapply(parms$aem$elements, function(i) ifelse(inherits(i, 'well'), reached_well(i, x = m[,1], y = m[,2]), FALSE), TRUE)
    rt <- as.numeric(!any(wls))
    return(rt)
  }
  # vectorized ODE
  # TODO check out other integration methods in {deSolve}
  get_paths <- Vectorize(function(x, y) {
    deSolve::lsoda(c(x, y), times = times, func = vxvy, parms = list(aem = aem, R = R),
                   events = list(root = TRUE),
                   rootfun = rootfun)
  }, SIMPLIFY = FALSE)

  # get paths and clean
  paths <- get_paths(x0, y0)
  paths.m <- lapply(paths, matrix, ncol = 3, dimnames = list(NULL, c('time', 'x', 'y')))
  if(length(paths) == 1) paths.m <- paths.m[[1]]
  class(paths.m) <- 'tracelines'
  return(paths.m)
}

#' Title
#'
#' @param tracelines
#'
#' @return
#' @export
#'
#' @examples
endpoints <- function(tracelines) {
  stopifnot('Supplied object is not of class {tracelines}' = inherits(tracelines, 'tracelines'))
  endp <- t(vapply(tracelines, function(i) i[nrow(i),], setNames(rep(0, 3), c('time', 'x', 'y'))))
  return(endp)
}


#' Title
#'
#' @param aem
#' @param well
#' @param time
#' @param npar
#' @param dt
#' @param as.poly
#' @param hull
#' @param ...
#' @param R
#'
#' @return
#' @export
#'
#' @examples
capzone <- function(aem, well, time, npar = 30, dt = time/100, as.poly = TRUE, hull = FALSE, R = 1, ...) {

  # define initial particle locations equally spaced at well screen circle
  alpha <- seq(0, 2*pi*(npar/(npar + 1)), length = npar)
  rw <- well$rw + 1e-12 # add small perturbation to prevent starting locations in root
  x <- rw*cos(alpha) + Re(well$zetaw)
  y <- rw*sin(alpha) + Im(well$zetaw)

  paths <- tracelines(aem, x0 = x, y0 = y, times = seq(0, time, dt), forward = FALSE, R = R, ...)
  if(as.poly) {
    endp <- endpoints(paths)[,-1]
    if(hull) endp <- endp[chull(endp),] # convex hull
    class(endp) <- 'capzone'
    return(endp)
  } else {
    return(paths)
  }
}
