
test_that("heads below aquifer base are dropped", {
  top <- 10; base <- 0
  uf <- uniformflow(100, 0.001, angle = 0)
  m <- aem(k = 10, top, base, n = 0.2, uf, type = 'variable')

  expect_warning(heads(m, 200, 0))
  expect_no_warning(heads(m, 200, 0, na.below = FALSE))
  expect_equal(heads(m, -200, 0), -heads(m, 200, 0, na.below = FALSE))

  base <- 1.2
  m <- aem(k = 10, top=4, base, n = 0.2, uf, type = 'variable')
  expect_equal(heads(m, -200, 0) - base, -(heads(m, 200, 0, na.below = FALSE) - base))

})
