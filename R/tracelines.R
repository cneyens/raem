

#' Check if point is on line
#'
#' [point_on_line()] checks if point p is on the line determined by points l0 and l1
#'
#' @param l0 numeric vector of length 3 with `x, y` and `z` coordinates at start of line
#' @param l1 numeric vector of length 3 with `x, y` and `z` coordinates at end of line
#' @param p numeric vector of length 3 with `x, y` and `z` coordinates of point to evaluate
#' @param width numeric, additional width of line. Defaults to 0 (no width).
#'
#' @return logical determining if point `p` is on line `l0-l1`
#' @noRd
#'
#' @examples
#'
#' p1 <- c(10, 20, 0)
#' p2 <- c(50, 40, 0)
#' m <- matrix(c(p1, p2), ncol = 3, byrow = TRUE)
#'
#' p3 <- c(30, 31, 0)
#'
#' plot(m, xlab = 'x', ylab = 'y')
#' lines(m)
#' points(matrix(p3, ncol = 3, byrow = TRUE), col = 'red')
#'
#'
#' point_on_line(p1, p2, p3)
#' point_on_line(p1, p2, c(40, 40, 0))
#' \dontrun{ # currently returns error
#' point_on_line(p1, p2, p3, width = 10)
#' }
#'
point_on_line <- function(l0, l1, p, width = 0, tol = 1e-6) {
  # line
  if(width == 0) {
    # https://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment?rq=3
    # crossproduct checks if new point is on the line as is extends to infinity
    # if yes, if dotproduct is positive and less than the squared length of the line, point falls on line between l0 and l1
    a <- l1 - l0
    b <- p - l1
    c <- p - l0
    crossp <- sum(c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] -
                      a[1] * b[3], a[1] * b[2] - a[2] * b[1])) # pracma::cross
    if(abs(crossp) > (width + tol)) return(FALSE)
    dotproduct <- sum(a * c) # pracma::dot
    if(dotproduct < 0) return(FALSE)
    sqL <- sum(a^2)
    if(dotproduct > sqL) return(FALSE)

    return(TRUE)
  } else {
    # basically rectangular polygon
    stop('width != 0 not yet supported', call. = FALSE)
  }
}

#' Check if particle has reached a line element
#'
#' @param line analytic element of class `linesink`
#' @param x x-coordinate of particle
#' @param y y-coordinate of particle
#' @param ... arguments passed to [point_on_line()]
#'
#' @details This assumes a fully penetrating line segment over the entire aquifer thickness.
#'
#' @return logical indicating if particle is on line
#' @noRd
#'
#' @examples
#' ls <- headlinesink(-75, 50, 100, 50, 10)
#' reached_line(ls, 50, 20)
#' reached_line(ls, 50, 50)
#'
reached_line <- function(line, x, y, ...) {
  l0 <- c(Re(line$z0), Im(line$z0), 0)
  l1 <- c(Re(line$z1), Im(line$z1), 0)
  p <- c(x, y, 0)
  width <- ifelse(is.null(line$width), 0, line$width) # TODO implement width in linesinks elements
  pol <- point_on_line(l0, l1, p, width = width, ...)
  return(pol)
}

#' Check if particle has reached inner annulus of well element
#'
#' @param well analytic element of class `well`
#' @param x x-coordinate of particle
#' @param y y-coordinate of particle
#' @param ... ignored
#'
#' @details This assumes a fully penetrating well over the entire aquifer thickness.
#'
#' @return logical indicating if particle has reached inner annulus of well
#' @noRd
#'
#' @examples
#' w <- well(50, 100, Q = 100, rw = 0.3)
#' reached_well(w, 51, 100)
#' reached_well(w, 50, 100.2)
#'
reached_well <- function(well, x, y, ...) {
  xw <- Re(well$zetaw)
  yw <- Im(well$zetaw)
  d <- sqrt((x - xw)^2 + (y - yw)^2) - well$rw
  # if(d <= 0) d <- 0
  return(d <= 0)
}

#' Check if point is vertically outside the saturated aquifer
#'
#' @param aem `aem` object
#' @param x x-coordinate of point
#' @param y y-coordinate of point
#' @param z z-coordinate of point
#' @param ... ignored
#'
#' @return A list with as first element a logical vector `outside` of length equal to `length(z)` indicating if the z-coordinate falls above
#'    the saturated groundwater level (i.e. the water-table for unconfined conditions or the aquifer
#'    top for confined conditions), or below the aquifer base. The second element `coords` is a vector with the z-coordinates reset to the saturated level
#'    or aquifer base if they were outside, or the original z-coordinates if no reset was necessary.
#' @noRd
#'
#' @examples
#' uf <- uniformflow(100, 0.001, 0)
#' rf <- constant(-1000, 0, 11)
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, uf, rf)
#' outside_vertical(m, x = c(-200, 0, 200), y = 0, z = c(12, -5, 10)) # confined
#' # TODO add unconfined example
outside_vertical <- function(aem, x, y, z, ...) {
  sat_lvl <- satthick(aem, x, y, as.grid = FALSE) + aem$base
  out <- z > sat_lvl | z < aem$base
  resetted <- ifelse(z > sat_lvl, sat_lvl, ifelse(z < aem$base, aem$base, z))
  return(list(outside = out, coords = resetted))
}

#' Compute tracelines of particles
#'
#' [tracelines()] tracks particle locations moving forward or backward with the advective groundwater flow
#' by numerically integrating the velocity vector. The resulting set of connected coordinates produce the
#' tracelines.
#'
#' @param aem `aem` object
#' @param x0 numeric vector, starting `x` locations of the particles
#' @param y0 numeric vector, starting `y` locations of the particles
#' @param z0 numeric vector, starting `z` locations of the particles
#' @param times numeric vector with the times at which locations should be registered
#' @param forward logical, should be forward (`TRUE`; default) or backward (`FALSE`) tracking be performed.
#' @param R numeric, retardation coefficient passed to [velocity()]. Defaults to 1 (no retardation).
#' @param tfunc function or list of functions with additional termination events for particles. See details. Defaults to `NULL`.
#' @param tol numeric tolerance used to define when particles have crossed a line element. Defaults to 0.1 length units.
#' @param ... ignored
#'
#' @details [deSolve::lsoda] is used to numerically integrate the velocity vector.
#'
#' Particles are terminated prematurely when they have reached the inner annulus of well elements, when they
#'    have crossed a line element or when they travel above the saturated aquifer level (i.e.
#'    the water-table for unconfined conditions or the aquifer top for confined conditions), or below the aquifer base.
#'    Note that these last two conditions can only occur in models with vertical flow components, i.e. models with area-sinks.
#'    The returned time value is the time of termination.
#'
#' The `tfunc` argument can be used to specify additional termination events. It is a function (or a list of functions) that
#'    takes arguments `t`, `coords` and `parms`. These are, respectively, a numeric value with the current tracking time,
#'    a numeric vector of length 3 with the current `x`, `y` and `z` coordinates of the particle, and a list with elements
#'    `aem` and `R`. It should return a single logical value indicating if the particle should terminate. See examples.
#'
#' If initial particle locations are above the saturated aquifer level, they are reset to this elevation with a warning.
#'    Initial particle locations below the aquifer base are reset at the aquifer base with a warning.
#'
#' Backward particle tracking is performed by reversing the flow field (i.e. multiplying the velocities with `-1`).
#'
#' @return [tracelines()] returns an object of class `tracelines` which is a list with length equal to the number of particles where each list element contains
#'    a matrix with columns `time`, `x`, `y` and `z` specifying the registered time and coordinates of the particle as is it tracked through the flow field.
#'
#' The final row represents either the location at the maximum `times` value or, if the particle terminated prematurely, the time and location of the termination.
#'
#' The matrices are ordered in increasing time. By connecting the coordinates, the tracelines can be produced.
#'
#' @rdname tracelines
#' @export
#' @seealso [capzone()]
#' @examples
#' k <- 10
#' top <- 10; base <- 0
#' n <- 0.2
#' R <- 5
#'
#' uf <- uniformflow(TR = 100, gradient = 0.001, angle = -10)
#' rf <- constant(TR, xc = -1000, yc = 0, hc = 10)
#'
#' m <- aem(k, top, base, n = n, uf, rf)
#'
#' x0 <- -200; y0 <- seq(-200, 200, 50)
#' times <- seq(0, 25 * 365, 365 / 20)
#' paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times)
#' endp <- endpoints(paths)
#'
#' xg <- seq(-500, 500, length = 100)
#' yg <- seq(-300, 300, length = 100)
#'
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(paths, add = TRUE, col = 'orange3')
#' points(endp[, c('x', 'y')])
#'
#' # Backward tracking
#' paths_back <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, R = R, forward = FALSE)
#' plot(paths_back, add = TRUE, col = 'forestgreen')
#'
#' # Termination at wells and linesinks
#' w1 <- well(200, 50, Q = 250)
#' w2 <- well(-200, -100, Q = 450)
#' ls <- headlinesink(-100, 100, 400, -300, 7)
#'
#' m <- aem(k, top, base, n = n, uf, rf, w1, w2, ls)
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(m, add = TRUE)
#'
#' x0 <- seq(-400, 400, 50); y0 <- 200
#' times <- seq(0, 5 * 365, 365 / 20)
#' paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times)
#' plot(paths, add = TRUE, col = 'orange3')
#'
#' # User-defined termination in rectangular zone
#' tzone <- cbind(x = c(-300, -200, -200, -300), y = c(150, 150, 100, 100))
#' termf <- function(t, coords, parms) {
#'   x <- coords[1]
#'   y <- coords[2]
#'   in_poly <- x <= max(tzone[,'x']) & x >= min(tzone[,'x']) &
#'              y <= max(tzone[,'y']) & y >= min(tzone[,'y'])
#'   return(in_poly)
#' }
#' paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, tfunc = termf)
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(m, add = TRUE)
#' polygon(tzone)
#' plot(paths, add = TRUE, col = 'orange3')
#'
#' # model with vertical flow due to area-sink
#' as <- areasink(xc = 0, yc = 0, N = 0.001, R = 1500)
#' m <- aem(k, top, base, n = n, uf, rf, w1, w2, as)
#'
#' # starting z0 locations are above aquifer top and will be reset to top with warning
#' paths <- tracelines(m, x0 = x0, y0 = y0, z0 = top + 0.5, times = times)
#'
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(m, add = TRUE)
#' plot(paths, add = TRUE, col = 'orange3')
#'
#' # plot vertical cross-section of traceline 4 along increasing y-axis (from south to north)
#' plot(paths[[4]][,c('y', 'z')], type = 'l')
#'
tracelines <- function(aem, x0, y0, z0, times, forward = TRUE, R = 1, tfunc = NULL, tol = 1e-1, ...) {

  if(!is.null(tfunc) & !is.list(tfunc)) tfunc <- list(tfunc)

  # reset initial locations if outside vertical domain
  outside_v_init <- outside_vertical(aem, x0, y0, z0)
  if(any(outside_v_init$outside)) warning('Resetting z0 values above saturated aquifer level or below aquifer base', call. = FALSE)
  z0 <- outside_v_init$coords

  # wrapper to obtain velocity
  vxvy <- function(t, coords, parms, ...) {
    m <- matrix(coords, ncol = 3, byrow = TRUE) # necessary ??
    v <- velocity(parms$aem, x=m[,1], y=m[,2], z=m[,3], R = parms$R, verbose = FALSE)
    v <- ifelse(is.na(v), 0, v) # Qz = NA if particle is above saturated aquifer or below aquifer base
    direction <- ifelse(forward, 1, -1)
    return(list(c(direction * v)))
  }

  # function to check for termination of particles
  rootfun <- function(t, coords, parms) {
    m <- matrix(coords, ncol = 3, byrow = TRUE) # necessary ??
    outside_v <- outside_vertical(parms$aem, m[,1], m[,2], m[,3])$outside
    wls <- vapply(parms$aem$elements, function(i) ifelse(inherits(i, 'well'), reached_well(i, x = m[,1], y = m[,2]), FALSE), TRUE)
    lls <- vapply(parms$aem$elements, function(i) ifelse(inherits(i, 'linesink'), reached_line(i, x = m[,1], y = m[,2], tol = tol), FALSE), TRUE)
    if(is.null(tfunc)) {
      tfn <- FALSE
    } else {
      tfn <- vapply(tfunc, function(i) do.call(i, list(t, coords, parms)), TRUE)
    }

    rt <- any(c(outside_v, wls, lls, tfn))
    return(as.numeric(!rt))
  }
  # vectorized ODE
  get_paths <- Vectorize(function(x, y, z) {
    deSolve::lsoda(c(x, y, z), times = times, func = vxvy, parms = list(aem = aem, R = R),
                   events = list(root = TRUE),
                   rootfun = rootfun)
  }, SIMPLIFY = FALSE)

  # get paths and clean
  paths <- get_paths(x0, y0, z0)
  paths.m <- lapply(paths, matrix, ncol = 4, dimnames = list(NULL, c('time', 'x', 'y', 'z')))
  # if(length(paths) == 1) paths.m <- paths.m[[1]] # return matrix instead of list of matrices when only one particle is tracked
  class(paths.m) <- 'tracelines'
  return(paths.m)
}

#' @description [endpoints()] obtains the final time and locations of tracked particles
#'
#' @param tracelines object of class `tracelines` as returned by [tracelines()]
#'
#' @return [endpoints()] returns a matrix with columns `time`, `x`, `y` and `z` specifying the final time and coordinates
#'     of the particles in the `tracelines` object.
#' @export
#' @rdname tracelines
endpoints <- function(tracelines, ...) {
  stopifnot('Supplied object is not of class {tracelines}' = inherits(tracelines, 'tracelines'))
  endp <- t(vapply(tracelines, function(i) i[nrow(i),], setNames(rep(0, 4), c('time', 'x', 'y', 'z'))))
  return(endp)
}


#' Calculate capture zone of a well element
#'
#' [capzone()] determines the capture zone of a well element in the flow field by performing backward
#' particle tracking until the requested time is reached.
#'
#' @param aem `aem` object
#' @param well well analytic element of class `well` or inherits from it.
#' @param time numeric, time of the capture zone
#' @param npar integer, number of particles to use in the backward tracking. Defaults to 30.
#' @param dt numeric, time discretization used in the particle tracking. Defaults `time / 100`.
#' @param zstart numeric vector with the starting elevations of the particles. Defaults to the base of the aquifer. See details.
#' @param as.poly logical, should the convex hull of the traces (default) be returned or the paths. See details.
#' @param ... additional arguments passed to [tracelines()].
#'
#' @details Backward particle tracking is performed using [tracelines()] and setting the `forward = FALSE`.
#'    Initial particle locations are computed by equally spacing `npar` locations at the well radius and repeating this
#'    for all `zstart` elevations. The total number of particles equals `npar * length(zstart)`.
#'
#' Note that different `zstart` values only have an effect in models with vertical flow components, i.e. models with area-sinks.
#'
#' If `as.poly = TRUE`, the convex hull of all particle locations is computed and returned. If `FALSE`,
#'   the output of [tracelines()] is returned.
#'
#' @return [capzone()] returns an object of class `capzone` if `as.poly = TRUE`, which is a matrix with columns `x` and `y`
#'    specifying the coordinates of the convex hull delineating the `time` capture zone of `well`.
#'    If `as.poly = FALSE`, the output of [tracelines()] is returned directly.
#' @rdname capzone
#' @export
#' @seealso [tracelines()]
#' @examples
#'
#' k <- 10
#' top <- 10; base <- 0
#' n <- 0.3
#'
#' uf <- uniformflow(TR = 100, gradient = 0.001, angle = -10)
#' rf <- constant(TR, xc = -1000, yc = 0, hc = 10)
#' w1 = well(200, 50, Q = 250)
#' w2 = well(-200, -100, Q = 450)
#'
#' m <- aem(k, top, base, n = n, uf, rf, w1, w2)
#'
#' cp5 <- capzone(m, w1, time = 5*365)
#' cp10 <- capzone(m, w2, time = 10*365, as.poly = FALSE)
#'
#' xg <- seq(-800, 800, length = 100)
#' yg <- seq(-500, 500, length = 100)
#'
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(cp5, add = TRUE)
#' plot(cp10, add = TRUE)
#'
#' # model with vertical flow components
#' as <- areasink(0, 0, N = 0.001, R = 1500)
#' m <- aem(k, top, base, n = n, uf, rf, w1, w2, as)
#'
#' # three different starting levels, returns convex hull of all particles
#' cp5 <- capzone(m, w1, time = 5*365, zstart = c(0, 5, 8))
#'
#' contour(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(cp5, add = TRUE)
#'
#' # one starting level at z = 8, smaller zone
#' cp5_z8 <- capzone(m, w1, time = 5*365, zstart = 8)
#' plot(cp5_z8, add = TRUE, col = 'orange3')
#'
capzone <- function(aem, well, time, npar = 30, dt = time / 100, zstart = aem$base, as.poly = TRUE, ...) {

  # define initial particle locations equally spaced at well screen circle
  # repeat for each elevation in zstart
  alpha <- seq(0, 2*pi*(npar/(npar + 1)), length = npar)
  rw <- well$rw + 1e-12 # add small perturbation to prevent starting locations in root
  x <- rw*cos(alpha) + Re(well$zetaw)
  y <- rw*sin(alpha) + Im(well$zetaw)
  c0 <- data.frame(x = rep(x, length(zstart)), y = rep(y, length(zstart)), z = rep(zstart, each = npar))

  paths <- tracelines(aem, x0 = c0$x, y0 = c0$y, z0 = c0$z, times = seq(0, time, dt), forward = FALSE, ...)
  if(as.poly) {
    pts <- do.call(rbind, paths)[,-1] # combine all particle locations in single matrix
    pts <- pts[chull(pts),] # convex hull
    class(pts) <- 'capzone'
    return(pts)
  } else {
    return(paths)
  }
}
