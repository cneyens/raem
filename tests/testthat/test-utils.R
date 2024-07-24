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
  expect_equal(setNames(c(0, -TR*i*l), c('positive', 'negative')),
               flow_through_line(m, x0, y0, x1, y1, split = TRUE))

  expect_equal(flow_through_line(m, x0, y0, x1, y1), -flow_through_line(m, x1, y1, x0, y0))


})

test_that('element_discharge works', {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  TR <- k * (top - base)

  rf <- constant(xc = -500, yc = 0, h = 20)
  uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
  w1 <- well(xw = 50, yw = 0, Q = 200)
  w2 <- well(xw = 0, yw = 100, Q = 400)
  hw <- headwell(xw = -100, yw = 0, hc = 7.5)
  hls <- headlinesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, hc = 8)
  as <- areasink(xc = 0, yc = 0, N = 0.0005, R = 500)
  m <- aem(k, top, base, n, rf, uf, w1, w2, hw, hls, as)

  expect_error(element_discharge(m)) # either name or type
  expect_error(element_discharge(m, name = 'hls', type = 'headlinesink')) # either name or type, not both
  expect_error(element_discharge(m, name = 'test')) # name not found
  expect_error(element_discharge(m, type = c('areasink', 'well'))) # only 1 type allowed
  expect_error(element_discharge(m, type = 'uniformflow')) # type should be allowed

  expect_equal(element_discharge(m, name = 'uf'), c('uf' = 0)) # return zero for unsupported types

  # area-sink
  expect_equal(element_discharge(m, type = 'areasink'), c('areasink' = -0.0005*pi*500^2))
  expect_equal(element_discharge(m, name = 'as'), c('as' = -0.0005*pi*500^2))

  # well
  expect_equal(element_discharge(m, type = 'well'), c('well' = w1$parameter + w2$parameter))
  expect_equal(element_discharge(m, 'w1'), c('w1' = w1$parameter))

  # line-sink
  expect_equal(element_discharge(m, type = 'headlinesink'), c('headlinesink' = m$elements$hls$parameter * hls$L))
  expect_equal(element_discharge(m, 'hls'), c('hls' = m$elements$hls$parameter * hls$L))

})
