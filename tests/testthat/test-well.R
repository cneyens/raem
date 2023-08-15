test_that('well potential and discharge is exact', {
  Qw <- 100
  w <- well(0, 0, Qw)

  # verify against potential form of Thiem eq. in polar coordinates
  # dis.thiem return x and y component
  pot.thiem <- function(x, y) Qw*log(sqrt(x^2 + y^2))/(2*pi)
  dis.thiem <- function(x, y) (-Qw/(2*pi*sqrt(x^2 + y^2)))*c(x/sqrt(x^2 + y^2), y/sqrt(x^2 + y^2))

  expect_equal(potential(w, 100, 0), pot.thiem(100, 0))
  expect_equal(potential(w, 0, 100), pot.thiem(0, 100))
  expect_equal(potential(w, 100, 100), pot.thiem(100, 100))

  expect_equal(c(Re(domega(w, 100, 0, 0)), -Im(domega(w, 100, 0, 0))), dis.thiem(100, 0))
  expect_equal(c(Re(domega(w, 0, 100, 0)), -Im(domega(w, 0, 100, 0))), dis.thiem(0, 100))
  expect_equal(c(Re(domega(w, 100, 100, 0)), -Im(domega(w, 100, 100, 0))), dis.thiem(100, 100))

})

test_that("well singularities are handled", {
  rw <- 0.3
  Qw <- 100
  w <- well(0, 0, Qw, rw = rw)

  # potential is the same everywhere within annulus
  # discharge depends on direction since it has x and y components
  expect_equal(potential(w, 0, 0), potential(w, rw, 0))
  expect_equal(potential(w, -0.1, 0.1), potential(w, rw, 0))

  expect_equal(domega(w, 0, 0, 0), domega(w, rw, 0, 0)) # Qy = 0

  # coordinates of -0.1, 0.1 projected onto annulus of radius rw
  x <- -0.1; y <- 0.1
  xrw <- rw*cos(atan2(y, x)); yrw <- rw*sin(atan2(y, x))
  expect_equal(domega(w, -0.1, 0.1, 0), domega(w, xrw, yrw, 0))

})

test_that("headwell works correctly", {

  k <- 10; top <- 10; base <- 0; n <- 0.2
  rf <- constant(-1000, 0, 10)
  hw <- headwell(xw = 200, yw = 100, hc = 8, rw = 0.3)
  m <- aem(k, top, base, n, rf, hw)

  expect_equal(heads(m, hw$xc, hw$yc), 8)

  hw2 <- headwell(xw = 200, yw = 100, hc = 8, rw = 0.3, xc = 100, yc = 100, rc = 0)
  m <- aem(k, top, base, n, rf, hw2)

  expect_equal(heads(m, 100, 100), 8)

  # resistance
  hc <- 8; res <- 20
  hw <- headwell(xw = 200, yw = 100, hc = hc, rw = 0.3, res = res)
  m <- aem(k, top, base, n, rf, hw, type = 'confined')
  Q <- (hc - heads(m, hw$xc, hw$yc))/(res) * (2*pi*hw$rw * satthick(m, hw$xc, hw$yc))
  expect_equal(Q, -m$elements$hw$parameter)

  # unconfined
  m <- aem(k, top, base, n, rf, hw, type = 'variable')
  Q <- (hc - heads(m, hw$xc, hw$yc))/(res) * (2*pi*hw$rw * satthick(m, hw$xc, hw$yc))
  expect_equal(Q, -m$elements$hw$parameter)
})
