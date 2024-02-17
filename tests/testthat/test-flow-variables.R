test_that("discharge works", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  i <- 0.001

  m <- aem(
    k = k, top = top, base = base, n = n,
    uf = uniformflow(TR = k * (top - base), gradient = i, angle = 0), type = "confined"
  )
  Q <- discharge(m, 50, 50, 0)
  expect_equal(Q[[1]], k * i * (top - base))
  expect_equal(Q[[2]], 0)
  expect_equal(colnames(Q), c("Qx", "Qy", "Qz"))

  N <- 0.0005
  z <- seq(10, 0, -0.5)
  m <- aem(
    k = k, top = top, base = base, n = n,
    uf = uniformflow(TR = k * (top - base), gradient = i, angle = 0),
    as = areasink(0, 0, N = N, R = 1000), type = "confined"
  )
  Q <- discharge(m, 50, 50, z)
  expect_equal(Q[, "Qz"], -N * (z - base))

  N <- 0.0005
  z <- seq(10, 0, -0.5)
  leakage <- -3 / (10 / 0.0001)
  m <- aem(
    k = k, top = top, base = base, n = n,
    uf = uniformflow(TR = k * (top - base), gradient = i, angle = 0),
    as = areasink(0, 0, N = N, R = 1000),
    asb = areasink(0, 0, N = leakage, R = 1000, loc = "base"), type = "confined"
  )
  Q <- discharge(m, 50, 50, z)
  sat <- satthick(m, 50, 50)
  expect_equal(Q[, "Qz"], (z - base) * (-N - leakage) + leakage * sat)
})

test_that("array names for flow are correct", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  i <- 0.001
  hc <- 20

  m <- aem(
    k = k, top = top, base = base, n = n,
    uf = uniformflow(TR = k * (top - base), gradient = i, angle = 0),
    rf = constant(-1000, 0, hc = hc)
  )
  x <- seq(-500, 500, length = 100)
  y <- seq(-100, 100, length = 100)
  z <- 2

  Q <- discharge(m, x, y, z)
  Qg <- discharge(m, x, y, z, as.grid = TRUE, magnitude = TRUE)
  expect_equal(colnames(Q), c("Qx", "Qy", "Qz"))
  expect_equal(dimnames(Qg)[[4]], c("Qx", "Qy", "Qz", "Q"))

  q <- darcy(m, x, y, z)
  qg <- darcy(m, x, y, z, as.grid = TRUE, magnitude = TRUE)
  expect_equal(colnames(q), c("qx", "qy", "qz"))
  expect_equal(dimnames(qg)[[4]], c("qx", "qy", "qz", "q"))

  v <- velocity(m, x, y, z)
  vg <- velocity(m, x, y, z, as.grid = TRUE, magnitude = TRUE)
  expect_equal(colnames(v), c("vx", "vy", "vz"))
  expect_equal(dimnames(vg)[[4]], c("vx", "vy", "vz", "v"))
})

test_that("discharge from phreatic curvature is correct", {
  top <- 10
  base <- 0
  k <- 10
  z <- 5

  grad <- 0.001
  TR <- k * (top - base)
  m <- aem(
    k = k, top = top, base = base, n = 0.2,
    uf = uniformflow(TR = TR, gradient = grad, angle = 0),
    rf = constant(-500, 0, 10), type = "variable"
  )
  h <- heads(m, 0, 0)
  Qx <- TR * grad # uniformflow: constant Q, independent of saturated thickness

  q <- darcy(m, 0, 0, z)
  # last term != grad because Q is constant in uniformflow; i.e. independent of satthick
  qz_exact <- (z / h) * (Qx / h) * (-Qx / (k * h))
  expect_equal(q[[1, "qz"]], qz_exact)

  # test with well ----
  Qw <- 100
  xw <- 250
  yw <- 0
  w <- well(xw, yw, Qw)
  m <- aem(k, top, base, n = 0.2, w, rf = constant(-500, 0, 8))

  x <- y <- 10
  z <- 5
  q <- darcy(m, x, y, z)
  h <- heads(m, x, y)

  dis.thiem <- function(x, y) {
    r <- sqrt((x - xw)^2 + (y - yw)^2)
    (-Qw / (2 * pi * r)) * c((x - xw) / r, (y - yw) / r)
  }
  Q_thiem <- dis.thiem(x, y)
  qz_exact <- (z / h) * ((Q_thiem[1] / h) * (-Q_thiem[1] / (k * h)) + (Q_thiem[2] / h) * (-Q_thiem[2] / (k * h)))

  expect_equal(q[[1, "qz"]], qz_exact)
})

test_that("discharge sets Qz to NA when z outside aquifer", {
  top <- 10
  base <- 0
  m <- aem(
    k = 10, top = top, base = base, n = 0.2,
    uf = uniformflow(TR = 100, gradient = 0.001, angle = 0)
  )
  expect_warning(Qt <- discharge(m, 0, 0, base - 1))
  expect_warning(Qb <- discharge(m, 0, 0, top + 1))
  expect_equal(Qt[[3]], NA_real_)
  expect_equal(Qb[[3]], NA_real_)
})

test_that("discharge returns array of correct dimensions", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  i <- 0.001
  hc <- 20

  m <- aem(
    k = k, top = top, base = base, n = n,
    uf = uniformflow(TR = k * (top - base), gradient = i, angle = 0),
    rf = constant(-1000, 0, hc = hc)
  )
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
  # expect_warning(expect_warning(expect_warning(expect_warning(discharge(m, x, y, z))))) # 6 warnings
  Q <- discharge(m, x, y, z, as.grid = TRUE)
  expect_equal(dim(Q), c(52, 100, 10, 3))
  Q <- discharge(m, x, y, z, as.grid = TRUE, magnitude = TRUE)
  expect_equal(dim(Q), c(52, 100, 10, 4))
})
