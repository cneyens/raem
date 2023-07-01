
constant <- function(TR, xc, yc, hc, ...) { # rename as reference ??
  cn <- element(0, 1)
  cn$xc <- xc
  cn$yc <- yc
  cn$pc <- TR*hc
  class(cn) <- c('constant', 'headequation', class(cn))
  return(cn)
}
omegainf.constant <- function(constant, x, y, ...) {
  omi <- as.complex(x*0 + 1)
  return(omi)
}

disc.constant <- function(constant, x, y, ...) {
  return(c(0 + 0*i))
}
