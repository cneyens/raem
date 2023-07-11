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

  expect_equal(c(discharge(w, 100, 0)), dis.thiem(100, 0))
  expect_equal(c(discharge(w, 0, 100)), dis.thiem(0, 100))
  expect_equal(c(discharge(w, 100, 100)), dis.thiem(100, 100))

})

test_that("well singularities are handled", {
  rw <- 0.3
  Qw <- 100
  w <- well(0, 0, Qw, rw = rw)

  # potential is the same everywhere within annulus
  # discharge depends on direction since it has x and y components
  expect_equal(potential(w, 0, 0), potential(w, rw, 0))
  expect_equal(potential(w, -0.1, 0.1), potential(w, rw, 0))

  expect_equal(discharge(w, 0, 0), discharge(w, rw, 0)) # Qy = 0

  # coordinates of -0.1, 0.1 projected onto annulus of radius rw
  x <- -0.1; y <- 0.1
  xrw <- rw*cos(atan2(y, x)); yrw <- rw*sin(atan2(y, x))
  expect_equal(discharge(w, -0.1, 0.1), discharge(w, xrw, yrw))

})
