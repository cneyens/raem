test_that('aem keeps names of element list', {

  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  TR <- k * (top - base)

  w <- well(xw = 50, yw = 0, Q = 200)
  rf <- constant(xc = -500, yc = 0, h = 20)
  uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)

  m <- aem(k, top, base, n, w, rf, uf)
  expect_equal(names(m$elements), c('w', 'rf', 'uf'))

  m <- aem(k, top, base, n, list(w, rf, uf))
  expect_equal(names(m$elements), NULL)

  m <- aem(k, top, base, n, list(well = w, constant = rf, flow = uf))
  expect_equal(names(m$elements), c('well', 'constant', 'flow'))

  m <- aem(k, top, base, n) |>
             add_element(rf, name = 'constant') |>
             add_element(w, name = 'well') |>
             add_element(uf, name = 'flow', solve = TRUE)
  expect_equal(names(m$elements), c('constant', 'well', 'flow'))
  expect_error(add_element(m, rf, name = 'constant'))
})
