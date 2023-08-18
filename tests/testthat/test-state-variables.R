
test_that('state-variables dimensions are correct and names are preserved', {
  top <- 10; base <- 0
  uf <- uniformflow(100, 0.001, angle = 0)
  rf <- constant(-500, 0, top*2)
  m <- aem(k = 10, top, base, n = 0.2, uf, rf)

  x <- seq(-20, 50, by = 10)
  expect_length(heads(m, x, y = 0), length(x))

  y <- seq(-10, 10, by = 10)
  expect_warning(heads(m, x, y))

  h <- heads(m, x, y, as.grid = TRUE)
  expect_equal(dim(h), c(length(y), length(x)))

  names(x) <- letters[seq_along(x)]
  expect_named(heads(m, x, 0), names(x))

  h <- heads(m, x, y, as.grid = TRUE)
  expect_named(h, NULL)

  y <- structure(x, names = LETTERS[seq_along(x)])
  h <- heads(m, x, y)
  expect_named(heads(m, x, y), names(x)) # names of x are preserved

  expect_named(heads(m, 0, y), names(y)) # names of x are preserved

})

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
