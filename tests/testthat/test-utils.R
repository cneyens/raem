test_that("saturated thickness is correct", {
  top <- 10
  base <- 0
  grad <- 0.001
  k <- 10
  TR <- k * (top-base)
  m <- aem(k = k, top = top, base = base, n = 0.2,
           uf = uniformflow(TR = TR, gradient = grad, angle = 0),
           rf = constant(-500, 0, 10.5), type = 'variable')
  xg <- seq(-500, 500, length = 100)
  h <- heads(m, xg, 0)
  sat <- satthick(m, xg, 0)
  expect_equal(sat[xg < 0], rep(top - base, sum(xg < 0)))
  expect_equal(sat[xg > 0], h[xg > 0])
})

test_that('directed flow works', {
  TR <- 100
  gradient <- 0.001
  angle <- -45
  k = 10; top = 10; base = 0; n = 0.2
  R <- 5
  rf <- constant(-1000, 0, hc = 10)
  uf <- uniformflow(TR = TR, gradient = gradient, angle = angle)
  m <- aem(k, top, base, n, rf, uf, type = 'confined')

  expect_equal(TR*gradient, dirflow(m, x = 40, y = -25, angle = angle))
  expect_equal(0, dirflow(m, x = 40, y = -25, angle = -angle))

  expect_equal(TR*gradient/(n*R*(top-base)), dirflow(m, x = 40, y = -25, angle = angle, flow = 'velocity', R = R))

  x <- 40:43
  y <- -(23:25)
  q <- dirflow(m, x = x, y = y, angle = angle, as.grid = TRUE, flow = 'darcy')
  expect_equal(dim(q), c(length(y), length(x)))


})

test_that('flow through line works', {
  TR <- 100; i <- 1e-3
  rf <- constant(-1000, 0, hc = 10)
  uf <- uniformflow(TR = TR, gradient = i, angle = -45)
  m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, uf)

  x0 <- -100
  y0 <- -100
  x1 <- 100
  y1 <- 100
  l <- sqrt((x1-x0)^2 + (y1-y0)^2)

  expect_equal(-TR*i*l, flow_through_line(m, x0, y0, x1, y1)) # negative Q is to the SE in this case
  expect_equal(flow_through_line(m, x0, y0, x1, y1), -flow_through_line(m, x1, y1, x0, y0))

})
