test_that('areasink complex discharge works', {
  TR <- 100; N <- 1e-3; R <- 700;
  m <- aem(TR, areasink(0, 0, N, R))
  exact <- function(r) ifelse(r <= R, N*r/2, N*R^2/(2*r))

  # exact at origin
  expect_equal(c(discharge(m, 0, 0)), c(exact(0), exact(0)))

  # exact signs
  expect_equal(c(discharge(m, 50, 0)), c(exact(50), exact(0)))
  expect_equal(c(discharge(m, -50, 0)), c(exact(-50), exact(0)))
  expect_equal(c(discharge(m, 0, 50)), c(exact(0), exact(50)))
  expect_equal(c(discharge(m, 0, -50)), c(exact(0), exact(-50)))

  # exact outside radius
  expect_equal(c(discharge(m, 0, 800)), c(exact(0), exact(800)))

})
