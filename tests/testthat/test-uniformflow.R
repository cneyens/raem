test_that("uniform flow has correct absolute values and angles", {
  TR <- 100
  gr <- 0.001
  uf <- uniformflow(TR, gradient = gr, angle = 0)
  x <- c(100, 200)
  y <- c(0, 0)
  h <- potential(uf, x, y) / TR
  expect_equal(diff(h), -gr * diff(x)) # flow from left to right, dh = i*dl

  angle <- -45
  uf <- uniformflow(TR, gradient = gr, angle = angle)
  x <- 50
  y <- 0
  r <- 100
  xn <- x + r * cos(angle * pi / 180)
  yn <- y + r * sin(angle * pi / 180)
  h <- potential(uf, c(x, y), c(xn, yn)) / TR
  expect_equal(diff(h), -gr * r)
})

test_that("uniform flow handles singularity", {
  uf <- uniformflow(1, 0.001, 0)
  expect_equal(omega(uf, 0, 0), c(0 + 0i))
  expect_equal(domega(uf, 0, 0), 0.001 + 0i)
})
