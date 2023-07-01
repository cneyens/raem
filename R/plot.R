


# plot
contour.aem <- function(aem, x, y, z = c('head', 'streamfunction', 'potential'), ...) {
  z <- match.arg(z)
  m <- switch(z,
              head = head(aem, x, y, as.grid = TRUE),
              streamfunction = streamfunction(aem, x, y, as.grid = TRUE),
              potential = potential(aem, x, y, as.grid = TRUE)
  )
  m <- t(m[dim(m)[1]:1,])
  return(contour(x, y, m, ...))
}
plot.element <- function(element, add = FALSE, pch = 16, cex = 0.75, ...) {
  if(any('well' %in% class(element))) {
    pts <- c(Re(element$zetaw), Im(element$zetaw))
    if(add) {
      return(points(pts[1], pts[2], pch = pch, cex = cex, ...))
    } else {
      return(plot(pts[1], pts[2], pch = pch, cex = cex, ...))
    }
  } else if(any(c('headlinesink', 'linesink') %in% class(element))) {
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
