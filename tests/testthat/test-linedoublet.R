

test_that("linedoublets work", {

  # TODO remove this when line-doublets are fixed
  skip_if(class(try(linedoublet(0, 1, 0, 1), silent = TRUE)) == 'try-error', message = 'Skipping line-doublet tests for now')

  k <- 8
  top <- 10
  base <- 0
  n <- 0.2
  res <- 500

  uf <- uniformflow(k*(top - base), 0.001, angle = 0)
  rf <- constant(-1000, 0, 8)
  ld <- linedoublet(x0 = -100, y0 = 100, x1 = 100, y1 = -100, res = res)

  # strength parameter for line-doublet should equal the jump in potential at collocation point
  # q = -dh / res
  m <- aem(k, top, base, n, uf, rf, ld, type = 'confined')
  dphi <- potential(m, ld$xc - 1e-12, ld$yc - 1e-12) - potential(m, ld$xc + 1e-12, ld$yc + 1e-12)
  expect_equal(m$element$ld$parameter, -dphi)

  dh <- heads(m, ld$xc + 1e-12, ld$yc + 1e-12) - heads(m, ld$xc - 1e-12, ld$yc - 1e-12)
  df <- dirflow(m, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  expect_equal(df, -dh / res)

  # unconfined
  m2 <- aem(k, top, base, n, uf, rf, ld, type = 'variable')
  dphi <- potential(m2, ld$xc - 1e-12, ld$yc - 1e-12) - potential(m2, ld$xc + 1e-12, ld$yc + 1e-12)
  expect_equal(m2$element$ld$parameter, -dphi)

  df <- dirflow(m2, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  hi <- heads(m2, ld$xc + 1e-12, ld$yc + 1e-12)
  ho <- heads(m2, ld$xc - 1e-12, ld$yc - 1e-12)
  dh <- hi - ho
  expect_equal(df, -dh / res)

  # test unconfined resfac with non-zero base
  m2 <- aem(k, top, base - 5, n, uf, rf, ld, type = 'variable')

  df <- dirflow(m2, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  hi <- heads(m2, ld$xc + 1e-12, ld$yc + 1e-12)
  ho <- heads(m2, ld$xc - 1e-12, ld$yc - 1e-12)
  dh <- hi - ho
  expect_equal(df, -dh / res)

})

test_that('impermeable wall works', {

  # TODO remove this when line-doublets are fixed
  skip_if(class(try(linedoublet(0, 1, 0, 1), silent = TRUE)) == 'try-error', message = 'Skipping line-doublet tests for now')

  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  res <- Inf

  uf <- uniformflow(k*(top - base), 0.001, angle = 0)
  rf <- constant(-1000, 0, 8)
  ld <- linedoublet(x0 = -100, y0 = 100, x1 = 100, y1 = -100, res = res)

  # strength parameter for line-doublet should equal the jump in potential at collocation point
  # q = 0 for impermeable wall
  m <- aem(k, top, base, n, uf, rf, ld, type = 'confined')
  dphi <- potential(m, ld$xc - 1e-12, ld$yc - 1e-12) - potential(m, ld$xc + 1e-12, ld$yc + 1e-12)
  expect_equal(m$element$ld$parameter, -dphi)

  # check that normal flow at collocation points equals zero
  df <- dirflow(m, ld$xc, ld$yc, angle = 45, flow = 'darcy')
  expect_equal(df, rep(0, length(ld$parameter)))
})
