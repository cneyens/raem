
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
  stopifnot(ndim <= 3)
  if(ndim == 2) {
    return(t(m[dims[1]:1,]))
  } else {
    return(aperm(m[dims[1]:1,,], c(2, 1, 3)))
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
  stopifnot(ndim <= 3)
  if(ndim == 2) {
    return(t(m[,dims[2]:1]))
  } else {
    return(aperm(m[,dims[2]:1,], c(2, 1, 3)))
  }
}

#' Display contours of a computed variable
#'
#' [contour.aem()] creates a contour plot of a computed variable, or adds the contour lines
#'     to and existing plot.
#'
#' @param aem `aem` object
#' @param x numeric x coordinates at which the values in `z` are computed These must be in ascending order.
#' @param y numeric x coordinates at which the values in `z` are computed These must be in ascending order.
#' @param z character indicating which variable to plot. Possible values are `head` (default),
#'    `streamfunction`, `potential`, `Q` (magnitude of discharge vector), `Qx`, `Qy`, `Qz`,
#'    (x, y and z components of the discharge vector) and `velocity`.
#' @param ... additional arguments passed to [contour()]
#'
#' @return Contour plot of variable `z`, or added contour plot of variable `z` to an
#'     existing plot.
#' @export
#'
#' @examples
#' TR <- 100
#' w <- well(xw = 50, yw = 0, Q = 200)
#' wi <- well(xw = 0, yw = 0, Q = -100)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
#' ml <- aem(TR = 100, w, wi, uf)
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#'
#' contour(ml, xg, yg, nlevels = 20, col = 'dodgerblue3')
#' contour(ml, xg, yg, 'streamfunction', nlevels = 20, col = 'orange3', drawlabels = FALSE, add = TRUE)
#'
contour.aem <- function(aem, x, y, z = c('head', 'streamfunction', 'potential'), ...) {
  z <- match.arg(z)
  m <- switch(z,
              head = head(aem, x, y, as.grid = TRUE),
              streamfunction = streamfunction(aem, x, y, as.grid = TRUE),
              potential = potential(aem, x, y, as.grid = TRUE)
  )
  m <- matrix_to_image(m)
  return(contour(x = x, y = y, z = m, ...))
}

#' Plot the location of an analytic element with a point or line geometry
#'
#' [plot.element()] plots the planar location of an analytic element with a point
#'    or line geometry or adds it to an existing plot.
#'
#' @param element analytic element of class `element` to plot. If not a point or line geometry, nothing is plotted.
#' @param add logical, should the plot be added to the existing plot? Defaults to `FALSE`.
#' @param pch numeric point symbol value, defaults to `16`.
#' @param cex numeric symbol size value, defaults to `0.75`.
#' @param ... additional arguments passed to [plot()]
#'
#' @return Plot of the planar location of the element, if any, or added to the existing plot.
#' @export
#'
#' @examples
#' ls <- linesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, sigma = 0.1)
#' w <- well(xw = 50, yw = 0, Q = 200)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#'
#' plot(ls)
#' plot(w, add = TRUE)
#' plot(uf) # empty
#'
plot.element <- function(element, add = FALSE, pch = 16, cex = 0.75, ...) {
  if(inherits(element, 'well')) {
    pts <- c(Re(element$zetaw), Im(element$zetaw))
    if(add) {
      return(points(pts[1], pts[2], pch = pch, cex = cex, ...))
    } else {
      return(plot(pts[1], pts[2], pch = pch, cex = cex, ...))
    }
  } else if(inherits(element, 'linesink')) {
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

#' Plot the locations of all elements in an `aem` object with point or line geometries
#'
#' [plot.aem()] plots the planar locations of all analytic elements with a point or line geometry
#'    in an `aem` object, or adds them to an existing plot.
#'
#' @param aem `aem` object
#' @param add logical, should the plot be added to the existing plot? Defaults to `FALSE`.
#' @param xlim numeric, plot limits along the x-axis. Required if `add = FALSE`.
#' @param ylim numeric, plot limits along the y-axis. Required if `add = FALSE`.
#' @param frame.plot logical, should a border be drawn around the plot. Defaults to `TRUE`.
#' @param ... additional arguments passed to [plot()].
#'
#' @return Plot of the planar location of all analytic elements, or added to the existing plot.
#' @export
#'
#' @examples
#' w <- well(xw = 50, yw = 0, Q = 200)
#' wi <- well(xw = 0, yw = 0, Q = -100)
#' uf <- uniformflow(gradient = 0.002, angle = -45, TR = 100)
#' ls <- linesink(-75, 50, 100, 50, sigma = 1)
#' ml <- aem(TR = 100, w, wi, uf, ls)
#'
#' plot(ml, xlim = c(-100, 100), ylim = c(-75, 75))
#'
#' xg <- seq(-100, 100, length = 500)
#' yg <- seq(-75, 75, length = 100)
#'
#' contour(ml, xg, yg, nlevels = 20, col = 'dodgerblue3')
#' plot(ml, add = TRUE)
#'
plot.aem <- function(aem, add = FALSE, xlim, ylim, frame.plot = TRUE, ...) {
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
