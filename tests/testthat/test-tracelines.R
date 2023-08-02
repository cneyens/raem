test_that("tracelines works for 2D flow", {

  k <- 10
  top <- 10; base <- 0
  n <- 0.2
  R <- 1.5
  i <- 0.001
  alpha <- -45
  TR <- k * (top - base)
  hc <- 15

  uf <- uniformflow(TR = TR, gradient = i, angle = alpha)
  rf <- constant(-1000, 0, hc = hc)
  m <- aem(k, top, base, n = n, uf, rf)

  x0 <- 0; y0 <- 0
  times <- seq(0, 10 * 365, 365 / 20)

  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, R = R)

  v <- -(k * i) / (n*R)
  dr <- v * times
  xp <- dr*sin(alpha * pi/180) + x0
  yp <- dr*cos(alpha * pi/180) + y0
  pts <- matrix(c(xp, yp), ncol = 2, dimnames = list(NULL, c('x', 'y')))
  expect_equal(paths[[1]][,c('x', 'y')], pts)

  # backward
  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, R = R, forward = FALSE)

  v <- -(k * i) / (n*R)
  dr <- -v * times
  xp <- dr*sin(alpha * pi/180) + x0
  yp <- dr*cos(alpha * pi/180) + y0
  pts <- matrix(c(xp, yp), ncol = 2, dimnames = list(NULL, c('x', 'y')))
  expect_equal(paths[[1]][,c('x', 'y')], pts)

  # multiple traces
  x0 <- 0; y0 <- seq(-100, 100, length = 10)
  times <- seq(0, 10 * 365, 365 / 20)

  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, R = R)

  v <- -(k * i) / (n*R)
  dr <- v * times
  l <- list()
  for(p in y0) {
    xp <- dr*sin(alpha * pi/180) + x0
    yp <- dr*cos(alpha * pi/180) + p
    pts <- matrix(c(times, xp, yp, rep(top, length(times))), ncol = 4, dimnames = list(NULL, c('time', 'x', 'y', 'z')))
    l[[length(l) + 1]] <- pts
  }
  class(l) <- 'tracelines'

  expect_equal(paths, l)

})

test_that('tracelines works for 3D flow', {

  k <- 10
  top <- 10; base <- 0
  b <- top - base
  n <- 0.2
  N <- 0.0005
  xas <- -1000
  R <- 2000

  as <- areasink(xas, 0, N = N, R = R)
  hls1 <- headlinesink(-500, -500, -500, 500, hc = 15)
  hls2 <- headlinesink(500, -500, 500, 500, hc = 10)
  rf <- constant(-1000, 0, 15)

  m <- aem(k, top, base, n = n, as, rf, hls1, hls2)

  x0 <- 0; y0 <- 0; times <- seq(0, 3*365, 365/20)
  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times)

  # Bakker & Post (2022), eq. 1.34, rearranged
  # vertical travel time in case of uniform recharge between two rivers
  zt <- function(t) b / exp(t*N/(n*b)) + base
  z_exact <- zt(times)

  expect_equal(paths[[1]][,'z'], z_exact)

})

test_that('termination in tracelines works', {

  k <- 10
  top <- 10; base <- 0
  n <- 0.2
  xw <- 200; yw <- 0; rw <- 0.3

  w <- well(xw, yw, 400, rw = rw)
  uf <- uniformflow(TR = 100, gradient = 0.001, angle = 0)
  rf <- constant(-1000, 0, 20)

  m <- aem(k, top, base, n = n, uf, w, type = 'confined')

  x0 <- -200; y0 <- seq(-100, 100, length = 10)
  times <- seq(0, 10 * 365, 365 / 20)

  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times)
  endp <- endpoints(paths)

  expect_true(all(abs(endp[,'x'] - xw) < 0.3))
  expect_true(all(abs(endp[,'y'] - yw) < 0.3))

  lsx <- xw
  lsy <- c(-500, 500)
  rf <- constant(-1000, 0, 15)
  ls <- headlinesink(lsx, lsy[1], lsx, lsy[2], hc = 12)
  m <- aem(k, top, base, n = n, uf, ls, rf, type = 'confined')

  tol <- 1e-1
  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, tol = tol)
  endp <- endpoints(paths)

  expect_equal(endp[,'x'], rep(lsx, nrow(endp)), tolerance = tol)

})

test_that('capzone works', {
  # Bakker & Post (2022), chapter 7.2, eq. 7.16
  n <- 0.25
  k <- 20
  top <- 15; base <- 0
  H <- top - base
  i <- -0.002
  U <- -k * H * i
  Q <- 400 # instead of 500
  rw <- 0.3

  trav_time <- function(x, y) {
    theta <- atan2(y, x)
    -n*H*x/U - (Q*n*H/(2*pi*U^2)) * (log(sin(theta - 2*pi*U*y/Q) / sin(theta)))
  }

  # aem
  w <- well(0, 0, Q, rw)
  uf <- uniformflow(k*H, gradient = -i, angle = 0)
  m <- aem(k, top, base, n, uf, w, type = 'confined')
  t <- 365*5
  cp5 <- capzone(m, w, time = t, as.poly = FALSE)
  endp <- endpoints(cp5)[-1,] # not include first value
  t_exact <- trav_time(endp[,'x'], endp[,'y'])

  expect_equal(t_exact, rep(t, nrow(endp)), tolerance = 0.1) # numeric tolerance

})
