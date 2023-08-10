
#' Convert a print style array to an image style array
#'
#' @param m 2D or 3D array in the conventional print layout style See details.
#' @details [image()] and related functions (such as [contour()]) require the x-axis to correspond to
#'     row number and the y-axis to the column number with column 1 at the bottom.
#'     This is a 90 degree contour-clockwise rotation of the conventional printed layout
#'     of a matrix, i.e. x-axis correspond to column numbers, y-axis correspond to row numbers,
#'     and value row 1 - column 1 is at the top-left.
#'
#' The adjustment is achieved by first reversing the 1st dimension (rows) and then transposing
#'     the rows and columns.
#'
#' @return The matrix `m` adjusted to conform to the style required by [image()] and [contour()]
#' @noRd
#' @seealso [matrix_to_image()], [image()]
matrix_to_image <- function(m) {
  dims <- dim(m)
  ndim <- length(dims)
  stopifnot(ndim <= 4)
  if(ndim == 2) {
    return(t(m[dims[1]:1,]))
  } else if(ndim == 3) {
    return(aperm(m[dims[1]:1,,, drop = FALSE], c(2, 1, 3)))
  } else {
    return(aperm(m[dims[1]:1,,,, drop = FALSE], c(2, 1, 3, 4)))
  }
}

#' Convert an image style array to a print style array
#'
#' @param m 2D or 3D array in the style required by [image()]. See details.
#' @details [image()] and related functions (such as [contour()]) require the x-axis to correspond to
#'     row number and the y-axis to the column number with column 1 at the bottom.
#'     This is a 90 degree contour-clockwise rotation of the conventional printed layout
#'     of a matrix, i.e. x-axis correspond to column numbers, y-axis correspond to row numbers,
#'     and value row 1 - column 1 is at the top-left.
#'
#' The adjustment is achieved by first reversing the 2nd dimension (columns) and then transposing
#'     the rows and columns.
#'
#' @return The matrix `m` adjusted to conform to the conventional printed layout.
#' @noRd
#' @seealso [matrix_to_image()], [image()]
image_to_matrix <- function(m) {
  dims <- dim(m)
  ndim <- length(dims)
  stopifnot(ndim <= 4)
  if(ndim == 2) {
    return(t(m[,dims[2]:1]))
  } else if(ndim == 3) {
    return(aperm(m[,dims[2]:1,, drop = FALSE], c(2, 1, 3)))
  } else {
    return(aperm(m[,dims[2]:1,,, drop = FALSE], c(2, 1, 3, 4)))
  }
}

#' Plot contours of a computed variable of the analytic element model
#'
#' @description [contours()] creates a contour plot of a variable computed by the analytic element
#'    model `aem`, or adds the contour lines to an existing plot.
#'
#' @param aem `aem` object
#' @param x numeric x coordinates at which the values in `z` are computed. These must be in ascending order.
#' @param y numeric x coordinates at which the values in `z` are computed. These must be in ascending order.
#' @param variable character indicating which variable to plot. Possible values are `heads` (default),
#'    `streamfunction` and `potential`.
#' @param asp the `y/x` aspect ratio, see [plot.window()]. Defaults to 1 (equal unit lengths).
#' @param ... additional arguments passed to [contour()].
#'
#' @details [contours()] is a wrapper around [contour()]. It obtains the values of `variable` at
#'    the supplied grid points and constructs the matrix supplied to [contour()] by reversing the rows and
#'    transposing the matrix (see also the documentation of [image()]).
#'
#' @export
#' @importFrom graphics contour
#' @seealso [aem()] [contour()] [image()] [heads()]
#' @examples
#'
#' w <- well(xw = 50, yw = 0, Q = 200)
#' wi <- well(xw = -200, yw = 0, Q = -100)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' rf <- constant(-1000, 0, hc = 10)
#' ml <- aem(k = 10, top = 10, base = 0, n = 0.2, w, wi, uf, rf)
#'
#' xg <- seq(-350, 200, length = 100)
#' yg <- seq(-125, 125, length = 100)
#'
#' contours(ml, xg, yg, nlevels = 20, col = 'dodgerblue3', labcex = 1)
#' contours(ml, xg, yg, 'streamfunction', nlevels = 20, col = 'orange3',
#'          drawlabels = FALSE, add = TRUE)
#'
#' # Not to be confused by contour()
#' try(
#' contour(ml, xg, yg, nlevels = 20, col = 'dodgerblue3', labcex = 1)
#' )
#'
#' # For image() or filled.contour()
#' library(graphics)
#' h <- heads(ml, xg, yg, as.grid = TRUE)
#' h_im <- t(h[dim(h)[1]:1,])
#' image(xg, yg, h_im, asp = 1)
#' contour(xg, yg, h_im, asp = 1, add = TRUE) # contours() is a wrapper for this
#' filled.contour(xg, yg, h_im, asp = 1)
#'
contours <- function(aem, x, y, variable = c('heads', 'streamfunction', 'potential'), asp = 1, ...) {
  variable <- match.arg(variable)
  if(!inherits(aem, 'aem')) stop('contours() should be called with an \'aem\' object', call. = FALSE)
  m <- switch(variable,
              heads = heads(aem, x, y, as.grid = TRUE),
              streamfunction = streamfunction(aem, x, y, as.grid = TRUE),
              potential = potential(aem, x, y, as.grid = TRUE)
  )
  m <- matrix_to_image(m)
  return(contour(x = x, y = y, z = m, asp = asp, ...))
}

#' @description [plot.element()] plots the location of an analytic element with point or line geometry.
#'
#' @param x `aem` object, or analytic element of class `element` to plot. If not a point or line geometry, nothing is plotted.
#' @param y ignored
#' @param add logical, should the plot be added to the existing plot? Defaults to `FALSE`.
#' @param pch numeric point symbol value, defaults to `16`.
#' @param cex numeric symbol size value, defaults to `0.75`.
#'
#' @details If the analytic element has a point geometry and has a collocation point
#'    (e.g. [headwell()]), that point is also plotted.
#'
#' A reference point (as created by [constant()]) is never plotted.
#'
#' @export
#' @importFrom graphics points lines plot.new
#' @rdname aem
#' @include aem.R
#' @examples
#' # Plotting ----
#' plot(ls)
#' plot(w, add = TRUE)
#' plot(uf) # empty
#'
plot.element <- function(x, y = NULL, add = FALSE, pch = 16, cex = 0.75, ...) {
  element <- x
  if(inherits(element, 'well')) {
    x <- Re(element$zetaw)
    y <- Im(element$zetaw)
    if(inherits(element, 'headwell')) {
      x[2] <- element$xc
      y[2] <- element$yc
    }
    if(add) {
      return(points(x, y, pch = pch, cex = cex, ...))
    } else {
      return(plot(x, y, pch = pch, cex = cex, ...))
    }
  } else if(inherits(element, 'linesink') || inherits(element, 'linedoublet')) {
    x <- c(Re(element$z0), Re(element$z1))
    y <- c(Im(element$z0), Im(element$z1))
    if(add) {
      return(lines(x, y, ...))
    } else {
      return(plot(x, y, type = 'l', ...))
    }
  } else {
    if(add) {
      invisible()
    } else {
      plot.new()
    }
  }
}

#'
#' @description [plot.aem()] plots the planar locations of all analytic elements with a point or line geometry
#'    in an `aem` object by calling [plot.element()] on them, or adds them to an existing plot.
#'
#' @param xlim numeric, plot limits along the x-axis. Required if `add = FALSE`.
#' @param ylim numeric, plot limits along the y-axis. Required if `add = FALSE`.
#' @param frame.plot logical, should a border be drawn around the plot. Defaults to `TRUE`.
#'
#' @export
#' @importFrom graphics frame plot.window axis box
#' @rdname aem
#' @include aem.R
#' @examples
#' plot(m, xlim = c(-500, 500), ylim = c(-250, 250))
#'
#' xg <- seq(-500, 500, length = 200)
#' yg <- seq(-250, 250, length = 100)
#'
#' contours(m, x = xg, y = yg, col = 'dodgerblue3', nlevels = 20)
#' plot(m, add = TRUE)
#'
plot.aem <- function(x, y = NULL, add = FALSE, xlim, ylim, frame.plot = TRUE, ...) {
  aem <- x
  if(add) {
    pl <- lapply(aem$elements, plot, add = add, ...)
    invisible(pl)
  } else {
    ln <- length(aem$elements)

    frame()
    plot.window(xlim, ylim)
    if(ln > 0) {
      for(i in seq_along(aem$elements)) plot(aem$elements[[i]], add = TRUE, ...)
    }
    axis(1)
    axis(2)
    box()
    invisible()
  }
}

#'
#' @param x object of class `tracelines`
#' @param y ignored
#' @param add logical, should the plot be added to the existing plot? Defaults to `FALSE`.
#' @param type character indicating what type of plot to draw. See [plot()]. Defaults to `'l'` (lines).
#' @param arrows logical indicating if arrows should be drawn using [arrows()]. Defaults to `FALSE`.
#' @param marker numeric, time interval at which to plot point markers. Defaults to `NULL` (no markers).
#' @param ... additional arguments passed to [plot()] or [arrows()].
#'
#' @export
#' @importFrom graphics arrows lines points
#' @rdname tracelines
#' @include tracelines.R
#' @examples
#'
#' # plot arrows
#' contours(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(paths, add = TRUE, col = 'orange3', arrows = TRUE, length = 0.05)
#'
#' # plot point markers every 2.5 years
#' contours(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
#' plot(paths, add = TRUE, col = 'orange3', marker = 2.5 * 365, pch = 20)
#'
plot.tracelines <- function(x, y = NULL, add = FALSE, type = 'l', arrows = FALSE, marker = NULL, ...) {
  if(!is.list(x)) {
    x <- list(x)
  }
  if(add) {
    if(arrows) {
      for(i in seq_along(x)) {
        arrows(x0 = x[[i]][1:(nrow(x[[i]])-1),2], y0 = x[[i]][1:(nrow(x[[i]])-1),3],
               x1 = x[[i]][2:nrow(x[[i]]),2], y1 = x[[i]][2:nrow(x[[i]]),3],
               ...)
      }
    } else {
      for(i in seq_along(x)) lines(x[[i]][,2], x[[i]][,3], type = type,...)
    }
  } else {
    if(arrows) {
      plot(x[[1]][,2], x[[1]][,3], type = 'l', ...)
      for(i in seq_along(x)) {
        arrows(x0 = x[[i]][1:(nrow(x[[i]])-1),2], y0 = x[[i]][1:(nrow(x[[i]])-1),3],
               x1 = x[[i]][2:nrow(x[[i]]),2], y1 = x[[i]][2:nrow(x[[i]]),3],
               ...)
      }
    } else {
      plot(x[[1]][,2], x[[1]][,3], type = type, ...)
      if(length(x) > 1) for(i in 2:length(x)) lines(x[[i]][,2], x[[i]][,3], type = type,...)
    }
  }

  # point markers
  if(!is.null(marker)) {
    if(is.numeric(marker)) {
      for(i in seq_along(x)) {
        ids <- which((x[[i]][, 'time'] %% marker) == 0)
        if(length(ids) > 0) points(x[[i]][ids, 'x'], x[[i]][ids, 'y'], ...)
      }
    } else {
      stop('marker should be numeric or NULL', call. = FALSE)
    }
  }

}

#'
#' @description [plot.capzone()] plots a `capzone` object as a polygon.
#'
#' @param x object of class `capzone`
#' @param y ignored
#' @param col colour of polygon. Defaults to a light grey with reduced opacity.
#' @param add logical, should the plot be added to the existing plot? Defaults to `FALSE`.
#' @param ... additional arguments passed to [tracelines()] for [capzone()] or additional arguments passed to [polygon()] when plotting.
#'
#' @export
#' @importFrom graphics polygon frame plot.window axis box
#' @rdname capzone
#' @include tracelines.R
#'
plot.capzone <- function(x, y = NULL, col = "#BEBEBE90", add = FALSE, ...) {
  if(add) {
    polygon(x, col = col, ...)
  } else {
    frame()
    plot.window(range(x[,1]), range(x[,2]))
    polygon(x, col = col, ...)
    axis(1); axis(2); box()
  }
}
