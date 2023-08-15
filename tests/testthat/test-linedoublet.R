test_that("linedoublet parameter is correct", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  res <- 500

  uf <- uniformflow(100, 0.001, angle = 0)
  rf <- constant(-1000, 0, 8)
  ld <- linedoublet(x0 = -100, y0 = 100, x1 = 100, y1 = -100, res = res)

  # strength parameter for line-doublet should equal the jump in potential at collocation point
  # q = dh / res
  m <- aem(k, top, base, n, uf, rf, ld, type = 'confined')
  dphi <- potential(m, ld$xc - 1e-12, ld$yc - 1e-12) - potential(m, ld$xc + 1e-12, ld$yc + 1e-12)
  expect_equal(m$element$ld$parameter, -dphi)

  dh <- heads(m, ld$xc + 1e-12, ld$yc + 1e-12) - heads(m, ld$xc - 1e-12, ld$yc - 1e-12)
  df <- dirflow(m, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  expect_equal(df, -dh * k / res) # ???
  # expect_equal(df, -dh / res)

  # unconfined
  skip("Linedoublets do not yet work properly for unconfined flow")
  m2 <- aem(k, top, base, n, uf, rf, ld, type = 'variable')
  dphi <- potential(m2, ld$xc - 1e-12, ld$yc - 1e-12) - potential(m2, ld$xc + 1e-12, ld$yc + 1e-12)
  expect_equal(m2$element$ld$parameter, -dphi)

  df <- dirflow(m2, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  dh <- heads(m2, ld$xc + 1e-12, ld$yc + 1e-12) - heads(m2, ld$xc - 1e-12, ld$yc - 1e-12)
  expect_equal(df, -dh * k / res) # ???
  # expect_equal(df, -dh / res)

})

# m1 <- m
# xg <- seq(-400, 400, length = 100); yg <- seq(-250, 250, length = 100)
# contours(m, xg, yg, col = 'dodgerblue3', nlevels=20)
# plot(m, add=T)
#
# contours(m2, xg, yg, col = 'dodgerblue3', nlevels=20)
# plot(m2, add=T)
#

