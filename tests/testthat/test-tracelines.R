test_that("tracelines works", {

  k <- 10
  top <- 10; base <- 0
  n <- 0.2
  R <- 1.5
  i <- 0.001
  alpha <- -45
  TR <- k * (top - base)
  uf <- uniformflow(TR = TR, gradient = i, angle = alpha)

  m <- aem(k, top, base, n = n, uf)

  x0 <- 0; y0 <- 0
  times <- seq(0, 10 * 365, 365 / 20)

  paths <- tracelines(m, x0 = x0, y0 = y0, z = top, times = times, R = R)

  v <- -(k * i) / (n*R)
  dr <- v * times
  xp <- dr*sin(alpha * pi/180) + x0
  yp <- dr*cos(alpha * pi/180) + y0
  pts <- matrix(c(xp, yp), ncol = 2, dimnames = list(NULL, c('x', 'y')))
  expect_equal(paths[[1]][,c('x', 'y')], pts)

})
