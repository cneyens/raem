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
