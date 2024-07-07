
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
#' contours(ml, xg, yg, nlevels = 20, col = 'dodgerblue', labcex = 1)
#' contours(ml, xg, yg, 'streamfunction', nlevels = 20, col = 'orange',
#'          drawlabels = FALSE, add = TRUE)
#'
#' # Not to be confused by contour()
#' try(
#' contour(ml, xg, yg, nlevels = 20, col = 'dodgerblue', labcex = 1)
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
#' @param pch numeric point symbol value, defaults to `16`. For a reference point, a value of 4 is used.
#' @param cex numeric symbol size value, defaults to `0.75`.
#' @param use.widths logical, if line elements with non-zero width are plotted, should they be plotted as polygons including
#'    the width (`TRUE`; default) or as infinitesimally thin lines (`FALSE`)?
#' @param col color of element. Defaults to `'black'`.
#'
#' @details If the analytic element has a point geometry and has a collocation point
#'    (e.g. [headwell()]), that point is also plotted with `pch = 1`.
#'
#' @export
#' @importFrom graphics points lines plot.new polygon frame plot.window axis box
#' @rdname aem
#' @include aem.R
#' @examples
#' # Plotting ----
#' plot(ls)
#' plot(w, add = TRUE)
#' plot(uf) # empty
#'
plot.element <- function(x,
                         y = NULL,
                         add = FALSE,
                         pch = 16,
                         cex = 0.75,
                         use.widths = TRUE,
                         col = 'black',
                         xlim,
                         ylim,
                         ...) {
  element <- x
  if(inherits(element, 'well')) {
    x <- Re(element$zetaw)
    y <- Im(element$zetaw)
    if(inherits(element, 'headwell')) {
      x[2] <- element$xc
      y[2] <- element$yc
      pch[2] <- 1
    }
    if(add) {
      return(points(x, y, pch = pch, cex = cex, col = col, ...))
    } else {
      return(plot(x, y, pch = pch, cex = cex, col = col, ...))
    }
  } else if(inherits(element, 'linesink') || inherits(element, 'linedoublet')) {
    x <- c(Re(element$z0), Re(element$z1))
    y <- c(Im(element$z0), Im(element$z1))
    width <- ifelse(is.null(element$width), 0, element$width)
    if(add) {
      if(use.widths && width > 0) {
        theta_norm <- atan2(Im(element$z1 - element$z0), Re(element$z1 - element$z0)) - pi/2
        xci <- x - 0.5 * width * cos(theta_norm)
        yci <- y - 0.5 * width * sin(theta_norm)
        xco <- x + 0.5 * width * cos(theta_norm)
        yco <- y + 0.5 * width * sin(theta_norm)
        poly <- matrix(c(xci[1], yci[1], xci[2], yci[2], xco[2], yco[2], xco[1], yco[1]),
                       byrow = TRUE, ncol = 2)
        return(polygon(poly, col = col, ...))
      } else {
        return(lines(x, y, col = col, ...))
      }
    } else {
      if(use.widths && width > 0) {
        theta_norm <- atan2(Im(element$z1 - element$z0), Re(element$z1 - element$z0)) - pi/2
        xci <- x - 0.5 * width * cos(theta_norm)
        yci <- y - 0.5 * width * sin(theta_norm)
        xco <- x + 0.5 * width * cos(theta_norm)
        yco <- y + 0.5 * width * sin(theta_norm)
        poly <- matrix(c(xci[1], yci[1], xci[2], yci[2], xco[2], yco[2], xco[1], yco[1]),
                       byrow = TRUE, ncol = 2)

        frame()
        plot.window(range(c(xci, xco)), range(c(yci, yco)))
        polygon(poly, col = col, ...)
        axis(1); axis(2); box()
      } else {
        return(plot(x, y, type = 'l', col = col, ...))
      }
    }
  } else if(inherits(element, 'areasink')) {
    alpha <- seq(0, 2 * pi, length = 100)
    if(add) {
      return(polygon(x = element$R * cos(alpha) + element$xc,
              y = element$R * sin(alpha) + element$yc,
              col = col,
              border = NA))
    } else {
      frame()
      plot.window(xlim, ylim)
      polygon(x = element$R * cos(alpha) + element$xc,
              y = element$R * sin(alpha) + element$yc,
              col = col,
              border = NA)
      axis(1)
      axis(2)
      box()
      invisible()
    }
  } else if(inherits(element, 'constant')) {
    if(add) {
      return(points(element$x, element$y, pch = 4, cex = cex, col = col, ...))
    } else {
      return(plot(element$x, element$y, pch = 4, cex = cex, col = col, ...))
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
#'
#' @details A reference point (as created by [constant()]) is never plotted when plotting the model
#'    as it is not a hydraulic feature. Area-sinks (as created by [areasink()] or [headareasink()])
#'    are also never plotted as they would clutter the plot. These elements can be plotted by
#'    calling `plot()` on them directly.
#'
#' @export
#' @importFrom graphics frame plot.window axis box
#' @rdname aem
#' @include aem.R
#' @seealso [contours()]
#' @examples
#' plot(m, xlim = c(-500, 500), ylim = c(-250, 250))
#'
#' xg <- seq(-500, 500, length = 200)
#' yg <- seq(-250, 250, length = 100)
#'
#' contours(m, x = xg, y = yg, col = 'dodgerblue', nlevels = 20)
#' plot(m, add = TRUE)
#'
plot.aem <- function(x, y = NULL, add = FALSE, xlim, ylim, ...) {
  aem <- x
  # drop areasinks and reference point
  el <- aem$elements[vapply(aem$elements, function(i) !(inherits(i, 'areasink') | inherits(i, 'constant')), TRUE)]
  if(add) {
    pl <- lapply(el, plot, add = add, ...)
    invisible(pl)
  } else {
    ln <- length(el)

    frame()
    plot.window(xlim, ylim)
    if(ln > 0) {
      for(i in seq_along(el)) plot(el[[i]], add = TRUE, ...)
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
#' @param marker numeric, time interval at which to plot point markers. Defaults to `NULL` (no markers). See details.
#' @param ... additional arguments passed to [plot()] or [arrows()].
#'
#' @details The `marker` value can be used to plot point markers at given time intervals, e.g. every 365 days (see examples).
#'    The x and y locations of each particle at the marked times are obtained by linearly interpolating from the computed particle locations.
#'
#' @export
#' @importFrom graphics arrows lines points
#' @importFrom stats approx
#' @rdname tracelines
#' @include tracelines.R
#' @examples
#'
#' # plot arrows
#' contours(m, xg, yg, col = 'dodgerblue', nlevels = 20)
#' plot(paths, add = TRUE, col = 'orange', arrows = TRUE, length = 0.05)
#'
#' # plot point markers every 2.5 years
#' contours(m, xg, yg, col = 'dodgerblue', nlevels = 20)
#' plot(paths, add = TRUE, col = 'orange', marker = 2.5 * 365, pch = 20)
#'
#' # plot point markers every 600 days
#' plot(paths, add = TRUE, col = 'forestgreen', marker = 600, pch = 1)
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
    if(is.numeric(marker) && length(marker) == 1) {
      for(i in seq_along(x)) {
        t <- seq(x[[i]][1, 'time'], x[[i]][nrow(x[[i]]), 'time'], by = marker)
        if(length(t) > 0) {
          xpts <- approx(x[[i]][, 'time'], x[[i]][, 'x'], xout = t)$y
          ypts <- approx(x[[i]][, 'time'], x[[i]][, 'y'], xout = t)$y
          points(xpts, ypts, ...)
        }
      }
    } else {
      stop('marker should be numeric of length 1 or NULL', call. = FALSE)
    }
  }

}
