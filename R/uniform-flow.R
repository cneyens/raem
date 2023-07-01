
uniformflow <- function(TR, gradient, angle, ...) {
  uf <- element(gradient * TR)
  uf$udir <- exp(-1i * (angle*pi/180))
  class(uf) <- c('uf', class(uf))
  return(uf)
}
omegainf.uf <- function(uf, x, y, ...) {
  omi <- -uf$udir * (x + y * 1i)
  return(omi)
}

disc.uf <- function(uf, x, y, ...) {
  om <- omega(uf, x, y, as.grid = FALSE, ...)
  zeta <- x + y*1i
  W <- -om/zeta
  return(W)
}
