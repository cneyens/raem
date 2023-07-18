
#' Compute the saturated thickness
#'
#' [satthick()] computes the saturated thickness of the aquifer of an `aem` object
#'     at the given x and y coordinates.
#'
#' @param aem `aem` object
#' @param x numeric x coordinates to evaluate at
#' @param y numeric y coordinates to evaluate at
#' @param as.grid logical, should a matrix of dimensions c(`length(y), length(x)`) be returned? Defaults to `FALSE`.
#' @param ... ignored
#'
#' @details If the aquifer is confined at `x` and `y`, the saturated thickness equals the aquifer thickness.
#'    For flow with variable saturated thickness, if the aquifer is unconfined at `x` and `y`, the saturated thickness
#'    it is calculated as the hydraulic head at `x` and `y` minus the aquifer base.
#'
#' @return [satthick()] returns a vector of `length(x)` (equal to `length(y)`) with the saturated thicknesses at `x` and `y`.
#'     If `as.grid = TRUE`, a matrix of dimensions `c(length(y), length(x))` described by
#'     marginal vectors `x` and `y` containing the saturated thicknesses at the grid points.
#' @export
#' @seealso [flow()], [state-variables()]
#' @examples
#' uf <- uniformflow(100, 0.001, 0)
#' rf <- constant(-1000, 0, 11)
#' m <- aem(k = 10, top = 10, base = 0, n = 0.2, uf, rf)
#'
#' satthick(m, x = c(-200, 0, 200), y = 0) # confined
#' satthick(m, x = seq(-500, 500, length = 100),
#'          y = seq(-250, 250, length = 100), as.grid = TRUE)
#' # TODO add unconfined example
#'
satthick <- function(aem, x, y, as.grid = FALSE, ...) {

  if(as.grid) {
    df <- expand.grid(x = x, y = y) # increasing x and y values
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }

  # TODO adjust for unconfined flow
  d <- aem$top - aem$base
  mb <- c(cbind(x = gx, y = gy, b = d)[,'b'], use.names = FALSE) # recycle x and y
  if(as.grid) {
    mb <- matrix(mb, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    mb <- image_to_matrix(mb)
  }
  return(mb)
}
