test_that("headlinesink works with resistance", {
  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  res <- 10
  width <- 3
  hc <- 6

  rf <- constant(-1000, 0, 10)
  hls <- headlinesink(x0 = -100, y0 = 100, x1 = 100, y1 = -100, hc = hc, res = res, width = width)
  m <- aem(k, top, base, n, rf, hls, type = "confined")

  # sigma = w * dh/c
  sigma <- width * (heads(m, 0, 0) - hc) / res
  expect_equal(sigma, m$elements$hls$parameter)

  # unconfined
  m <- aem(k, top, base, n, rf, hls, type = "variable")
  sigma <- width * (heads(m, 0, 0) - hc) / res
  expect_equal(sigma, m$elements$hls$parameter)

  # test unconfined resfac with non-zero base
  m <- aem(k, top, base - 5, n, rf, hls, type = "variable")
  sigma <- width * (heads(m, 0, 0) - hc) / res
  expect_equal(sigma, m$elements$hls$parameter)
})
