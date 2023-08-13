test_that("linedoublet parameter is correct", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2

  uf <- uniformflow(100, 0.001, angle = 0)
  rf <- constant(-1000, 0, 10)
  ld <- linedoublet(x0 = -100, y0 = 100, x1 = 100, y1 = -100, res = 500) # s = -16.63472582

  # strength parameter for line-doublet should equal the jump in potential at collocation point
  m <- aem(k, top, base, n, uf, rf, ld, type = 'confined')
  dphi <- potential(m, ld$xc + 1e-12, ld$yc + 1e-12) - potential(m, ld$xc - 1e-12, ld$yc - 1e-12)
  expect_equal(m$element$ld$parameter, dphi)
})
