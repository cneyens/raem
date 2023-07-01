
linesink <- function(x0 = 1, y0 = 0, x1 = 1, y1 = 1, sigma = 1, ...) {
  ls <- element(sigma)
  ls$z0 <- x0 + y0 * 1i
  ls$z1 <- x1 + y1 * 1i
  ls$L <- abs(ls$z1 - ls$z0)
  class(ls) <- c('linesink', class(ls))
  return(ls)
}
omegainf.linesink <- function(linesink, x, y, ...) {
  zeta <- x + y * 1i
  Z <- (2 * zeta - (linesink$z0 + linesink$z1)) / (linesink$z1 - linesink$z0)
  tol <- 1e-12
  zp1 <- ifelse(abs(Z + 1) < tol, tol, Z + 1)
  zm1 <- ifelse(abs(Z - 1) < tol, tol, Z - 1)
  omi <- linesink$L / (4*pi) * (zp1*log(zp1) - zm1*log(zm1))
  return(omi)
}
headlinesink <- function(TR, x0 = 0, y0 = 0, x1 = 1, y1 = 1, hc = 1, ...) {
  hls <- linesink(x0 = x0, y0 = y0, x1 = x1, y1 = y1, sigma = 0)
  hls$xc <- 0.5*(x0 + x1)
  hls$yc <- 0.5*(y0 + y1)
  hls$pc <- TR*hc
  hls$nunknowns <- 1
  class(hls) <- c('headlinesink', 'headequation', class(hls))
  return(hls)
}
