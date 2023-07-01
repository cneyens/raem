
areasink <- function(xc = 0, yc = 0, N = 0.001, R = 100, ...) {
  as <- element(N)
  as$xc <- xc
  as$yc <- yc
  as$R <- R
  class(as) <- c('areasink', class(as))
  return(as)
}
omegainf.areasink <- function(areasink, x, y, ...) {
  Rs <- areasink$R
  r <- sqrt((x - areasink$xc)^2 + (y - areasink$yc)^2)
  phi <- ifelse(r < Rs, -0.25*(r^2 - Rs^2), -Rs^2/2 * log(r/Rs)) + 0i
  return(phi)
}

disc.areasink <- function(areasink, x, y, ...) { # TODO, absolute values not correct
  # Haitjema 1995, eq. 5.30, for exfiltration
  Rs <- areasink$R
  r <- sqrt((x - areasink$xc)^2 + (y - areasink$yc)^2)
  # Bakker & Post, 2022, eq. 6.39 & 6.40
  Qr <- areasink$parameter * ifelse(r <= Rs, r/2, Rs^2/(2*r))
  W <- (Qr*(x - areasink$xc)/(2*pi*r^2)) + (-Qr*(y - areasink$yc)/(2*pi*r^2))*1i # not sure about this
  W[r < 1e-15] <- 0 + 0i # continuous across boundary
  return(W)
}
