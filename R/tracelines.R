
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
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
tracelines <- function(aem, x0, y0, times, forward = TRUE, ...) {

  # wrapper to obtain velocity
  vxvy <- function(t, y, parms, ...) {
    m <- matrix(y, ncol = 2, byrow = TRUE) # necessary ??
    v <- velocity(parms$aem, x=m[,1], y=m[,2], n=parms$n, b=parms$b)
    return(list(c(v)))
  }

  # TODO other functions when particles should terminate
  rootfun <- function(t, y, parms) {
    m <- matrix(y, ncol = 2, byrow = TRUE) # necessary ??
    wls <- vapply(parms$aem$elements, function(i) ifelse(inherits(i, 'well'), reached_well(i, x = m[,1], y = m[,2]), FALSE), TRUE)
    rt <- as.numeric(!any(wls))
    return(rt)
  }
  # vectorized ODE
  get_paths <- Vectorize(function(x, y) {
    deSolve::lsoda(c(x, y), times = times, func = vxvy, parms = list(aem = aem, n = 0.2, b = 10),
                   events = list(root = TRUE),
                   rootfun = rootfun)
  }, SIMPLIFY = FALSE)

  # get paths and clean
  paths <- get_paths(x0, y0)
  paths.m <- lapply(paths, matrix, ncol = 3, dimnames = list(NULL, c('time', 'x', 'y')))
  if(length(paths) == 1) paths.m <- paths.m[[1]]
  class(paths.m) <- c('tracelines', class(paths.m))
  return(paths.m)
}
