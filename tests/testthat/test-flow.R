test_that("discharge works", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  i <- 0.001

  m <- aem(k = k, top = top, base = base, n = n,
           uf = uniformflow(TR = k*(top-base), gradient = i, angle = 0))
  Q <- discharge(m, 50, 50, 0)
  expect_equal(Q[[1]], k*i*(top-base))
  expect_equal(Q[[2]], 0)
  expect_equal(colnames(Q), c('Qx', 'Qy', 'Qz'))

  N <- 0.0005
  z <- seq(10, 0, -0.5)
  m <- aem(k = k, top = top, base = base, n = n,
           uf = uniformflow(TR = k*(top-base), gradient = i, angle = 0),
           as = areasink(0, 0, N = N, R = 1000))
  Q <- discharge(m, 50, 50, z)
  expect_equal(Q[,'Qz'], -N*(z - base))
})

test_that('discharge sets Qz to NA when z outside aquifer', {
  top <- 10
  base <- 0
  m <- aem(k = 10, top = top, base = base, n = 0.2,
           uf = uniformflow(TR = 100, gradient = 0.001, angle = 0))
  expect_warning(Qt <- discharge(m, 0, 0, base - 1))
  expect_warning(Qb <- discharge(m, 0, 0, top + 1))
  expect_equal(Qt[[3]], NA_real_)
  expect_equal(Qb[[3]], NA_real_)
})

test_that('discharge returns array of correct dimensions', {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  i <- 0.001

  m <- aem(k = k, top = top, base = base, n = n,
           uf = uniformflow(TR = k*(top-base), gradient = i, angle = 0))
  x <- 50
  y <- 50
  z <- c(0, 10)
  Q <- discharge(m, x, y, z)
  expect_equal(dim(Q), c(2, 3))
  Q <- discharge(m, x, y, z, magnitude = TRUE)
  expect_equal(dim(Q), c(2, 4))

  x <- seq(-500, 500, length = 100)
  y <- seq(-100, 100, length = 100)
  z <- 10
  Q <- discharge(m, x, y, z)
  expect_equal(dim(Q), c(100, 3))

  y <- seq(-100, 100, length = 52)
  z <- seq(10, 0, length = 10)
  expect_warning(expect_warning(expect_warning(discharge(m, x, y, z)))) # 3 warnings
  Q <- discharge(m, x, y, z, as.grid = TRUE)
  expect_equal(dim(Q), c(52, 100, 10, 3))
  Q <- discharge(m, x, y, z, as.grid = TRUE, magnitude = TRUE)
  expect_equal(dim(Q), c(52, 100, 10, 4))

})
