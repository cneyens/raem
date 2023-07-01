
equation <- function(element, aem, ...) {
  if(!inherits(element, 'headequation')) stop('element should be of class headequation', call. = FALSE)
  row <- vector(mode = 'numeric')
  rhs <- element$pc
  xc <- element$xc
  yc <- element$yc
  for(i in aem$elements) {
    if(i$nunknowns == 1) {
      row[length(row)+1] <- potinf(i, xc, yc)
    } else {
      rhs <- rhs - potential(i, xc, yc)
    }
  }
  return(list(row, rhs))
}
