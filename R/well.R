
well <- function(xw = 0, yw = 0, Q = 1, rw = 0.3, ...) {
  well <- element(Q)
  well$zetaw <- xw + 1i * yw
  well$rw <- rw
  class(well) <- c('well', class(well))
  return(well)
}
omegainf.well <- function(well, x, y, ...) {
  zminzw <- x + 1i * y - well$zetaw
  zminzw <- ifelse(abs(zminzw) < well$rw, well$rw, zminzw)
  omi <- 1/(2*pi) * log(zminzw)
  return(omi)
}

disc.well <- function(well, x, y, ...) {
  zeta <- x + y*1i
  W <- -well$parameter/(2*pi * (zeta - well$zetaw))
  return(W)
}


headwell <- function(TR, xw = 0, yw = 0, hw = 0, rw = 0.3, ...) {
  hwe <- well(xw = xw, yw = yw, Q = 0, rw = rw)
  hwe$xc <- xw + rw
  hwe$yc <- yw
  hwe$pc <- TR * hw
  hwe$nunknowns <- 1
  class(hwe) <- c('headwell', 'headequation', class(hwe))
  return(hwe)
}
