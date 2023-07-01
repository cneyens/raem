
aem <- function(TR, ...) {
  l <- list(...)
  # if(length(l) == 1 && inherits(l[[1]], 'list')) l <- l[[1]]
  # https://stackoverflow.com/questions/72894519/deparse-and-substitute-on-ellipsis-to-get-names-of-parameters
  names(l) <- sapply(substitute(...()), deparse)
  if(any(vapply(l, function(i) !inherits(i, 'element'), FALSE))) stop('All supplied objects should be analytic elements', call. = FALSE)
  aem <- list(TR = TR, elements = l, solved = FALSE)
  class(aem) <- 'aem'

  if(any(vapply(aem$elements, function(i) inherits(i, 'headequation'), TRUE))) {
    aem <- solve(aem)
  }

  return(aem)
}

omega <- function(...) UseMethod('omega')
potential <- function(...) UseMethod('potential')
streamfunction <- function(...) UseMethod('streamfunction')
potinf <- function(...) UseMethod('potinf')
omegainf <- function(...) UseMethod('omegainf')

omega.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y) # increasing x and y values
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  om <- lapply(aem$elements, omega, x = gx, y = gy)
  om <- colSums(do.call(rbind, om))
  # om <- 0 + 0i
  # for(i in aem$elements) om <- om + omega(i, gx, gy)

  if(as.grid) {
    om <- matrix(om, nrow = length(x), ncol = length(y)) # as used by {image} or {contour}. NROW and NCOL are switched
    om <- t(om[,dim(om)[2]:1])
  }
  return(om)
}
potential.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  pt <- Re(omega(aem, x, y, as.grid = as.grid))
  return(pt)
}
streamfunction.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  sf <- Im(omega(aem, x, y, as.grid = as.grid))
  return(sf)
}
head.aem <- function(aem, x, y, as.grid = FALSE, ...) { # rename as heads ???
  hd <- potential(aem, x, y, as.grid = as.grid) / aem$TR
  return(hd)
}


solve.aem <- function(aem, ...) {
  if(!any(vapply(aem$elements, function(i) inherits(i, 'headequation'), TRUE))) {
    stop('Model should contain at least 1 headequation element in order to be solved', call. = FALSE)
  }

  nun <- vapply(aem$elements, function(i) i$nunknowns, 1)
  esolve <- aem$elements[which(nun == 1)]
  nunknowns <- length(esolve)
  m <- matrix(0, nrow = nunknowns, ncol = nunknowns)
  rhs <- rep(0, nunknowns)
  for(irow in 1:nunknowns) {
    eq <- equation(esolve[[irow]], aem)
    m[irow,] <- eq[[1]]
    rhs[irow] <- eq[[2]]
  }
  solution <- solve(m, rhs)
  for(irow in 1:nunknowns) esolve[[irow]]$parameter <- solution[irow]
  aem$elements[which(nun == 1)] <- esolve
  aem$solved <- TRUE
  return(aem)
}

omega.element <- function(element, x, y, ...) {
  om <- element$parameter * omegainf(element, x, y)
  return(om)
}

potential.element <- function(element, x, y, ...) {
  pt <- Re(omega(element, x, y))
  return(pt)
}

potinf.element <- function(element, x, y, ...) {
  pti <- Re(omegainf(element, x, y))
  return(pti)
}
element <- function(p, un = 0, ...) {
  el <- list()
  el$parameter <- p
  el$nunknowns <- un
  class(el) <- c('element')
  return(el)
}
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
headwell <- function(TR, xw = 0, yw = 0, hw = 0, rw = 0.3, ...) {
  hwe <- well(xw = xw, yw = yw, Q = 0, rw = rw)
  hwe$xc <- xw + rw
  hwe$yc <- yw
  hwe$pc <- TR * hw
  hwe$nunknowns <- 1
  class(hwe) <- c('headwell', 'headequation', class(hwe))
  return(hwe)
}
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
add_element <- function(aem, element, name = NULL, solve = FALSE, ...) {
  if(!inherits(aem, 'aem')) stop('{aem} should be of class aem', call. = FALSE)
  if(!inherits(element, 'element')) stop('{element} should be of class element', call. = FALSE)
  if(is.null(aem$elements)) aem$elements <- list()
  aem$elements[[length(aem$elements) + 1]] <- element
  if(!is.null(name)) {
    if(name %in% names(aem$elements)) stop('element ', '\'', name, '\'', ' already exists', call. = FALSE)
    names(aem$elements)[length(aem$elements)] <- name
  }
  if(inherits(element, 'headequation') & solve) aem <- solve(aem)
  return(aem)
}
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

# complex discharge W at given point
# superposition of complex discharges of all elements

# complex discharge functions
disc <- function(...) UseMethod('disc')
# disc.element # not necessary

disc.uf <- function(uf, x, y, ...) {
  om <- omega(uf, x, y, as.grid = FALSE, ...)
  zeta <- x + y*1i
  W <- -om/zeta
  return(W)
}
disc.well <- function(well, x, y, ...) {
  zeta <- x + y*1i
  W <- -well$parameter/(2*pi * (zeta - well$zetaw))
  return(W)
}
disc.constant <- function(constant, x, y, ...) {
  return(c(0 + 0*i))
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
disc.aem <- function(aem, x, y, as.grid = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  w <- lapply(aem$elements, disc, x = gx, y = gy)
  w <- colSums(do.call(rbind, w))
  # w <- 0 + 0i
  # for(i in aem$elements) w <- w + disc(i, gx, gy)

  if(as.grid) {
    w <- matrix(w, nrow = length(x), ncol = length(y))  # as used by {image} or {contour}. NROW and NCOL are switched
    w <- t(w[,dim(w)[2]:1])
  }
  return(w)
}

# discharge vector
# Qx = Re(W); Qy = -Im(W)
# returns matrix with two columns (Qx, Qy) when as.grid = FALSE
# otherwise 3D matrix with x,y,2 dimensions ([,,1] = Qx, [,,2] = Qy)
# no {discharge} functions per element TYPE (e.g. 'uf' or 'well'): they dispatch on element class
# magnitude = TRUE returns an additional column or layer with the magnitude of the vector computed as
# Q = sqrt(Qx^2 + Qy^2)
disvec <- function(...) UseMethod('disvec')
disvec.element <- function(element, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  if(as.grid) {
    df <- expand.grid(x = x, y = y)
    gx <- df$x
    gy <- df$y
  } else {
    gx <- x
    gy <- y
  }
  W <- disc(element, x, y, ...)
  Qx <- Re(W)
  Qy <- -Im(W)
  if(magnitude) {
    qv <- c(Qx, Qy, sqrt(Qx^2 + Qy^2))
    ndim <- 3
    nms <- c('Qx', 'Qy', 'Q')
  } else {
    qv <- c(Qx, Qy)
    ndim <- 2
    nms <- c('Qx', 'Qy')
  }

  if(as.grid) {
    Q <- array(qv, dim = c(length(x), length(y), ndim), dimnames = list(NULL, NULL, nms)) # as used by {image} or {contour}. NROW and NCOL are switched
    Q <- aperm(Q[,dim(Q)[2]:1,], c(2,1,3))
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}
disvec.aem <- function(aem, x, y, as.grid = FALSE, magnitude = FALSE, ...) {
  W <- disc(aem, x, y, as.grid = as.grid, ...)
  Qx <- Re(W)
  Qy <- -Im(W)
  if(magnitude) {
    qv <- c(Qx, Qy, sqrt(Qx^2 + Qy^2))
    ndim <- 3
    nms <- c('Qx', 'Qy', 'Q')
  } else {
    qv <- c(Qx, Qy)
    ndim <- 2
    nms <- c('Qx', 'Qy')
  }
  if(as.grid) {
    Q <- array(qv, dim = c(dim(Qx), ndim), dimnames = list(NULL, NULL, nms)) # as used by {image} or {contour}. NROW and NCOL are switched
    Q <- aperm(Q[,dim(Q)[2]:1,], c(2,1,3))
  } else {
    Q <- matrix(qv, ncol = ndim, dimnames = list(NULL, nms))
  }
  return(Q)
}


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

# examples
TR = 100
w <- well(xw = 50, yw = 0, Q = 200)
wi <- well(xw = 0, yw = 0, Q = -100)
uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
ml <- aem(TR = 100, w, wi, uf)

xg <- seq(-100, 100, length = 500)
yg <- seq(-75, 75, length = 100)

contour(ml, xg, yg, 'potential', nlevels = 20, col = 'dodgerblue3')
contour(ml, xg, yg, 'stream', nlevels = 20, col = 'orange3', drawlabels = FALSE, add = TRUE)
plot(ml, add = TRUE)

ml <- aem(TR = 100)
w1 <- well(xw = 0, yw = 0, Q = 100, rw = 0.3)
w2 <- headwell(xw = 400, yw = 300, rw = 0.3, hw = 20, TR = 100)
uf <- uniformflow(TR = 100, gradient = 0.002, angle = 0)
rf <- constant(TR = 100, xc = 400, yc = 0, hc = 22)
ml <- aem(TR = 100, w1, w2, uf, rf)

xg <- seq(-200, 600, length = 100)
yg <- seq(-100, 500, length = 100)
contour(ml, xg, yg, levels = seq(20, 40, 0.2), col = 'dodgerblue3', xlab = 'x (m)', ylab = 'y (m)', grid = grid())
contour(ml, xg, yg, 'stream', levels = seq(-200, 100, ml$elements$w2$parameter/20), col = 'orange3', add = TRUE, drawlabels = FALSE)
plot(ml, add = TRUE)

ls1 <- linesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, sigma = 0.1)
ml <- aem(TR = 100, ls1)

xg <- seq(-400, 400, length = 100)
yg <- seq(-300, 300, length = 100)
h <- head(ml, xg, yg, as.grid = TRUE)
xs <- seq(-200, 200, length = 100)
ys <- seq(-150, 150, length = 100)
hs <- head(ml, xs, ys)

contour(ml, xg, yg, col = 'dodgerblue3', xlab = 'x (m)', ylab = 'y (m)', grid = grid())
contour(ml, xg, yg, 'stream', col = 'orange3', add = TRUE, drawlabels = FALSE)
plot(ml, add = TRUE)
plot(sqrt((xs+200)^2 + (ys+150)^2), hs, type = 'l')

TR <- 100
hriver <- 10
ml <- aem(TR = TR)
ml <- add_element(ml, constant(TR = TR, xc = 0, yc = 1000, hc = hriver + 2), 'rf')
ml <- add_element(ml, headwell(TR = TR, xw = 0, yw = 100, rw = 0.3, hw = hriver - 2), 'headwell', solve = FALSE)
xls <- seq(-1600, 1600, length = 101)
xls <- c(seq(-1600, -400, 200), seq(-350, 400, 50), seq(600, 1601, 200))
yls <- 50*sin(pi*xls/400)
for(i in seq_len(length(xls)-1)) {
  hls <- headlinesink(TR = TR, x0 = xls[[i]], y0 = yls[[i]], x1 = xls[[i+1]], y1 = yls[[i + 1]], hc = hriver)
  ml <- add_element(ml, hls, paste("headlinesink", i, sep = '_'), solve = FALSE)
}

ml <- solve(ml)

xg1 <- seq(-1800, 1800, length = 100)
yg1 <- seq(-1200, 1200, length = 101)
xg2 <- seq(-400, 400, length = 100)
yg2 <- seq(-100, 400, length = 100)
h1 <- head(ml, xg1, yg1, as.grid = TRUE)
h2 <- head(ml, xg2, yg2, as.grid = TRUE)

contour(ml, xg1, yg1, nlevels = 10, col = 'dodgerblue3')
plot(ml, add = TRUE)
contour(ml, xg2, yg2, nlevels = 20, col = 'dodgerblue3')
plot(ml, add = TRUE)
contour(ml, xg2, yg2, z = 'stream', col = 'orange3', add = TRUE)

TR = 100
rf <- constant(TR = TR, xc = 0, yc = 0, hc = 20)
as1 <- areasink(xc = -500, yc = 0, N = 0.001, R=500)
as2 <- areasink(xc = 500, yc = 0, N = -0.001, R = 500)
ml <- aem(TR = TR, rf, as1, as2)
xg <- seq(-1500, 1500, length = 100)
yg <- seq(-800, 800, length = 101)
h <- head(ml, xg, yg, as.grid = T)

contour(ml, xg, yg, nlevels = 10, col = 'dodgerblue3')
filled.contour(xg, yg, t(h[dim(h)[1]:1,]), color.palette = function(n) hcl.colors(n, 'viridis'))
plot(xg, h[,51], type = 'l')

Q <- disvec(ml, xg, y = 100)
plot(xg, Q[,1], type = 'l')
lines(xg, Q[,2], type = 'l', col = 'red')

hds <- head(ml, c(-450, -450.00001), y = 0)
grad <- diff(hds)/0.00001
Qx <- disvec(ml, -450, y = 0)[,1]

xls0 <- c(0, 100, 200, 400, 600, 800, 1000, 1100, 1200)
yls0 <- c(200, 200, 100, 100, 0, 0, 100, 300, 450)
hls0 <- seq(39, 40.4, length = 8)
xls1 <- c(0, 0, 200, 400, 600, 800, 1000, 1100, 1200)
yls1 <- c(200, 400, 600, 600, 700, 700, 750, 800, 850)
hls1 <- seq(39, 40.4, length = 8)

TR <- 100
ml <- aem(TR) |>
  add_element(constant(TR, xc = 0, yc = 800, hc = 39.5), 'rf') |>
  add_element(well(xw = 500, yw = 250, Q = 100), 'w0') |>
  add_element(well(xw = 800, yw = 500, Q = 100), 'w1') |>
  add_element(areasink(xc = 600, yc= 400, N = 0.001, R = 700))
for(i in seq_along(hls0)) {
  ml <- add_element(ml, headlinesink(TR, x0 = xls0[i], y0 = yls0[i], x1 = xls0[i+1], y1 = yls0[i+1], hc = hls0[i]),
                    name = paste('hls0', i, sep = '_'))
}
for(i in seq_along(hls1)) {
  ml <- add_element(ml, headlinesink(TR, x0 = xls1[i], y0 = yls1[i], x1 = xls1[i+1], y1 = yls1[i+1], hc = hls1[i]),
                    name = paste('hls1', i, sep = '_'))
}
ml <- solve(ml)

xg <- seq(-100, 1300, length = 100)
yg <- seq(-100, 900, length = 100)
h <- head(ml, xg, yg, as.grid = TRUE)

contour(ml, xg, yg, levels = seq(38, 41, 0.1), asp = 1, col = 'dodgerblue3', panel.first = grid())
plot(ml, add = TRUE)

filled.contour(xg, yg, h, color.palette = function(n) hcl.colors(n, 'viridis'),
               plot.axes = {plot(ml, add = TRUE)})

plot(ml, xlim = range(xg), ylim = range(yg), panel.first = grid(), asp = 1)

h <- head(ml, xg, y = 400)
plot(xg, h, type = 'l')

TR = 100
b = 10 # thickness
w <- well(xw = 100, yw = 150, Q = 200)
uf <- uniformflow(gradient = 0.002, angle = 50, TR = TR)
rf <- constant(TR, -100, -1000, 10)
ml <- aem(TR = 100, uf, rf, w)

xg <- seq(0, 300, length = 100)
yg <- seq(0, 300, length = 100)
contour(ml, xg, yg)

Q <- disvec(ml, xg, yg, as.grid = TRUE, magnitude = TRUE)
contour(xg, yg, Q[,,3], nlevels = 100)
