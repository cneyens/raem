test_that("discharge works", {
  m <- aem(k = 10, top = 10, base = 0, n = 0.2,
           uf = uniformflow(TR = 100, gradient = 0.001, angle = 0))
  Q <- discharge(m, 50, 50, 0)
  expect_equal(Q[[1]], 0.1)
  expect_equal(Q[[2]], 0)
  expect_equal(colnames(Q), c('Qx', 'Qy', 'Qz'))
})
