
# note that these tests only check that there are no errors or warnings
# produced while plotting. They don't actually verify if the plot is correct
# It is recommended to use expect_snapshot() for that, or use the {vdiffr} package.
# Both approaches compare the created plot output to a previously created, known output file

test_that("contour plots work", {
  k <- 10
  top <- 10
  base <- -10
  n <- 0.2

  uf <- uniformflow(TR = k * (top - base), gradient = 0.002, angle = 0)
  rf <- constant(x = -1000, y = 0, h = 8)
  w <- well(x = 0, y = 0, Q = 400)

  m <- aem(k, top, base, n, uf, rf, w)

  xg <- seq(-500, 500, length = 100)
  yg <- seq(-300, 300, length = 100)

  expect_invisible(contours(m, xg, yg, nlevels = 20, col = 'dodgerblue'))
  expect_invisible(contours(m, xg, yg, variable = 'streamfunction', nlevels = 20, col = 'orange', drawlabels = TRUE))
  expect_invisible(contours(m, xg, yg, 'potential', levels = seq(1000, 1400, 50), col = 'forestgreen'))

  expect_error(contours(list(x = 'A'), xg, yg))

})

test_that("plotting elements work", {
  k <- 10
  top <- 10
  base <- -10
  n <- 0.2

  uf <- uniformflow(TR = k * (top - base), gradient = 0.002, angle = 0)
  rf <- constant(x = -550, y = 0, h = 8)
  w <- well(x = 0, y = 0, Q = 400)
  hw <- headwell(xw = 200, yw = 0, hc = 4)
  as <- areasink(x = 0, y = 0, R = 500, N = 0.2 / 365)
  hls <- headlinesink(x0 = -300, y0 = -200, x1 = -300, y1 = 200, h = 6.8)

  m <- aem(k, top, base, n, uf, rf, w, hw, as, hls)

  xg <- seq(-500, 500, length = 100)
  yg <- seq(-300, 300, length = 100)

  expect_invisible(contours(m, xg, yg, nlevels = 20, col = 'dodgerblue'))
  expect_invisible(plot(m, add = TRUE))
  expect_invisible(plot(as, add = TRUE, col = adjustcolor('grey60', alpha = 0.6)))
  expect_invisible(plot(rf, add = TRUE))

  hls <- headlinesink(x0 = -300, y0 = -200, x1 = -300, y1 = 200, h = 6.8, width = 10)

  m <- aem(k, top, base, n, uf, rf, w, hw, as, hls)
  expect_invisible(plot(m, xlim = c(-600, 600), ylim = c(-300, 300)))

  expect_invisible(plot(w, xlim = c(-600, 600), ylim = c(-300, 300)))
  expect_invisible(plot(hls, xlim = c(-600, 600), ylim = c(-300, 300)))
  expect_invisible(plot(hls, xlim = c(-600, 600), ylim = c(-300, 300), use.widths = FALSE))
  expect_invisible(plot(as, xlim = c(-600, 600), ylim = c(-300, 300), col = adjustcolor('grey60', alpha = 0.6)))
  expect_invisible(plot(uf)) # empty
  expect_invisible(plot(rf))

})

test_that('plotting tracelines works', {

  k <- 10
  top <- 10
  base <- -10
  n <- 0.2

  uf <- uniformflow(TR = k * (top - base), gradient = 0.002, angle = 0)
  rf <- constant(x = -550, y = 0, h = 8)

  m <- aem(k, top, base, n, uf, rf)
  paths <- tracelines(m, x0 = -400, y0 = seq(-200, 200, length = 3), z0 = base, times = seq(0, 5*365, 365/10))

  xg <- seq(-500, 500, length = 100)
  yg <- seq(-300, 300, length = 100)

  expect_invisible(contours(m, xg, yg, nlevels = 20, col = 'dodgerblue'))
  expect_invisible(plot(paths, add = TRUE))
  expect_invisible(plot(paths, add = TRUE, marker = 163, col = 'forestgreen'))
  expect_error(plot(paths, add = TRUE, marker = c(10, 20), col = 'forestgreen'))
  expect_invisible(plot(paths, add = TRUE, arrows = TRUE, col = 'orange'))

  expect_invisible(plot(paths, add = FALSE, xlim = c(-500, 500), ylim = c(-300, 300)))
  expect_invisible(plot(paths, add = FALSE, xlim = c(-500, 500), ylim = c(-300, 300), arrows = TRUE))

})
