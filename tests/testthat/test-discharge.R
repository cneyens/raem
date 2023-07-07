test_that("discharge works", {
  m <- aem(TR = 100, uf = uniformflow(TR = 100, gradient = 0.001, angle = 0))
  Q <- discharge(m, 50, 50)
  expect_equal(Q[[1]], 0.1)
  expect_equal(Q[[2]], 0)
  expect_equal(colnames(Q), c('Qx', 'Qy'))
})
