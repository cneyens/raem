test_that('areasink complex discharge works', {
  N <- 1e-3; R <- 700; xc <- 600; yc <- 400
  m <- aem(k = 10, top = 10, base = 0, n = 0.2,
           areasink(xc, yc, N, R))
  exact <- function(r) ifelse(r <= R, N*r/2, N*R^2/(2*r))

  # exact at origin
  expect_equal(domega(m, xc, yc), c(exact(0) - exact(0)*1i))

  # exact signs
  expect_equal(domega(m, xc + 50, yc + 0), c(exact(50) - exact(0)*1i))
  expect_equal(domega(m, xc - 50, yc + 0), c(exact(-50) - exact(0)*1i))
  expect_equal(domega(m, xc, yc + 50), c(exact(0) - exact(50)*1i))
  expect_equal(domega(m, xc, yc - 50), c(exact(0) - exact(-50)*1i))

  # exact outside radius
  expect_equal(domega(m, xc, yc + 800), c(exact(0) - exact(800)*1i))

})

test_that('headareasink works', {

  rf <- constant(-1000, 0, 15)
  has <- headareasink(0, 0, hc = 20, R = 1500)
  m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, has, type = 'confined')
  as <- areasink(0, 0, N = m$elements$has$parameter, R=1500)
  m2 <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, as, type = 'confined')

  expect_equal(heads(m, seq(0, 100, length = 11), 0), heads(m2, seq(0, 100, length = 11), 0))

  hc <- 20
  res <- 800
  has <- headareasink(0, 0, hc = hc, R = 1500, res = res)
  m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, has, type = 'confined')

  qz <- (hc - heads(m, has$xc, has$yc))/res
  expect_equal(m$elements$has$parameter, qz)

  hc <- 14
  res <- 800
  has <- headareasink(0, 0, hc = hc, R = 1500, res = res)
  m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, has, type = 'confined')

  qz <- (hc - heads(m, has$xc, has$yc))/res
  expect_equal(m$elements$has$parameter, qz)

  # unconfined
  hc <- 9
  res <- 800
  rf <- constant(-1000, 0, 8)
  has <- headareasink(0, 0, hc = hc, R = 1500, res = res)
  m <- aem(k = 10, top = 10, base = 0, n = 0.2, rf, has, maxiter = 100)

  qz <- (hc - heads(m, has$xc, has$yc)) / res
  expect_equal(m$elements$has$parameter, qz)

})
