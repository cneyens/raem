test_that('areasink complex discharge works', {
  N <- 1e-3; R <- 700;
  m <- aem(k = 10, top = 10, base = 0, n = 0.2,
           areasink(0, 0, N, R))
  exact <- function(r) ifelse(r <= R, N*r/2, N*R^2/(2*r))

  # exact at origin
  expect_equal(domega(m, 0, 0), c(exact(0) - exact(0)*1i))

  # exact signs
  expect_equal(domega(m, 50, 0), c(exact(50) - exact(0)*1i))
  expect_equal(domega(m, -50, 0), c(exact(-50) - exact(0)*1i))
  expect_equal(domega(m, 0, 50), c(exact(0) - exact(50)*1i))
  expect_equal(domega(m, 0, -50), c(exact(0) - exact(-50)*1i))

  # exact outside radius
  expect_equal(domega(m, 0, 800), c(exact(0) - exact(800)*1i))

})
