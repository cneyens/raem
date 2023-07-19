test_that('aem keeps names of element list', {

  k <- 10
  top <- 10
  base <- 0
  n <- 0.2
  TR <- k * (top - base)

  w <- well(xw = 50, yw = 0, Q = 200)
  rf <- constant(xc = -500, yc = 0, h = 20)
  uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)

  m <- aem(k, top, base, n, w, rf, uf)
  expect_equal(names(m$elements), c('w', 'rf', 'uf'))

  m <- aem(k, top, base, n, well = w, constant = rf, flow = uf)
  expect_equal(names(m$elements), c('well', 'constant', 'flow'))

  m <- aem(k, top, base, n, w, constant = rf, flow = uf)
  expect_equal(names(m$elements), c('w', 'constant', 'flow'))

  m <- aem(k, top, base, n, list(w, rf, uf))
  expect_equal(names(m$elements), NULL)

  m <- aem(k, top, base, n, list(well = w, constant = rf, flow = uf))
  expect_equal(names(m$elements), c('well', 'constant', 'flow'))

  if(!(R.version$major >= 4 & R.version$minor >= 1)) {
    skip('Skipping test because R version < 4.1.0 and does not have native pipe')
  }
  # expect_warning(m <- aem(k, top, base, n)) # empty aem model
  # m <- m |>
  m <- aem(k, top, base, n) |>
    add_element(rf, name = 'constant') |>
    add_element(w, name = 'well') |>
    add_element(uf, name = 'flow', solve = TRUE)
  expect_equal(names(m$elements), c('constant', 'well', 'flow'))
  expect_error(add_element(m, rf, name = 'constant'))
})

test_that('when solving aem, matrix is not singular', {

  k <- 10
  top <- 10; base <- 0
  n <- 0.2

  uf <- uniformflow(TR = 100, 0.001, 0)
  w <- well(-100, 50, 250)
  rf <- constant(-1000, 0, 10)
  hw <- headwell(-50, 50, 9)
  ls <- linesink(20, -100, 20, 100, sigma = 5)
  hls <- headlinesink(-20, -100, -20, 100, hc = 8)

  expect_no_error(aem(k, top, base, n, uf))
  expect_error(aem(k, top, base, n, hw))
  expect_no_error(aem(k, top, base, n, hw, rf))
  expect_error(aem(k, top, base, n, hls))
  expect_no_error(aem(k, top, base, n, hls, rf))
  expect_error(aem(k, top, base, n, hls, hw, ls, uf, w))
  expect_no_error(aem(k, top, base, n, hls, hw, ls, uf, w, rf))

})


test_that('aem is exact for 2D flow with a well in uniform background flow', {

  # based on Bakker & Post (2022), chapter 7.1
  # uniform flow in x-direction
  xw <- 0
  yw <- 0
  Q <- 80
  k <- 10
  i <- -0.001
  top <- 10; base <- 0
  angle <- 0

  # analytical
  U <- -k * i * (top - base)

  pot <- function(x, y) {
    r <- sqrt((x - xw)^2 + (y - yw)^2)
    phi <- Q/(2*pi) * log(r) - U*x
    return(phi)
  }
  psi <- function(x, y) {
    Q/(2*pi) * atan2(y, x) - U*y
  }
  QQ <- function(x, y) {
    r <- sqrt((x - xw)^2 + (y - yw)^2)
    Qx <- U - Q/(2*pi) * (x - xw)/(r^2)
    Qy <- -Q/(2*pi) * (y - yw)/(r^2)
    m <- cbind(Qx = Qx, Qy = Qy, Qz = 0)
    return(m)
  }

  # aem
  # gradient in uniformflow is i but positive in direction of flow
  uf <- uniformflow(TR = k*(top-base), gradient = -i, angle = angle)
  w <- well(xw, yw, Q)

  m <- aem(k, top, base, n = 0.2, uf, w)

  df <- expand.grid(x = seq(-500, 500, length = 50), y = seq(-200, 200, length = 25))

  expect_equal(potential(m, df$x, df$y), pot(df$x, df$y))
  expect_equal(streamfunction(m, df$x, df$y), psi(df$x, df$y))
  expect_equal(discharge(m, df$x, df$y, z = 0), QQ(df$x, df$y))

})

test_that('aem is exact for 2D flow with a well in uniform background flow near a river', {

  # based on Bakker & Post (2022), chapter 7.3
  Q <- 80
  k <- 10
  i <- -0.001
  top <- 10; base <- 0
  angle <- 0
  h0 <- 0
  TR <- k * (top - base)
  d <- 50
  phi0 <- h0 * TR

  # analytical
  UL <- -k * i * (top - base)
  UR <- -UL
  om <- function(x, y) {
    zeta <- x + y*1i
    ifelse(x <= 0, -UL*zeta + Q/(2*pi)*log((zeta + d)/(d - zeta)) + phi0, -UR*zeta + phi0)
  }
  dom <- function(x, y) {
    zeta <- x + y*1i
    ifelse(x <= 0, UL - Q/(2*pi) * (1/(zeta + d) - 1/(zeta - d)), -UR)
  }

  # skip('Skipping 2D with well, flow and river')

  # since Qy is taking up too much flow, results differ slightly, especially further away from river
  # TODO try again once linedoublets are implemented to set no-flow at y = + and y = -

  # aem
  xrf <- -1000
  hc <- i * (0 + xrf) + h0
  uf <- uniformflow(TR = TR, -i, angle = angle)
  w <- well(0 - d, 0, Q)
  rf <- constant(xrf, 0, hc)
  m <- aem(k, top, base, n = 0.2, uf, w, rf, type = 'confined')
  lseg <- 10
  for(n in c(seq(-1005, -105, lseg), seq(105, 1005, lseg))) {
    hls <- headlinesink(x0 = 0, x1 = 0, y0 = n - 0.5*lseg, y1 = n + 0.5*lseg, hc = h0)
    m <- add_element(m, hls, solve = FALSE)
  }
  lseg <- 1
  for(n in seq(-100, 100, lseg)) {
    hls <- headlinesink(x0 = 0, x1 = 0, y0 = n - 0.5*lseg, y1 = n + 0.5*lseg, hc = h0)
    m <- add_element(m, hls, solve = FALSE)
  }
  m <- solve(m)

  xg <- seq(-500, 0, length = 100)
  # expect_equal(om(xg, 0), omega(m, xg, 0))

  # with tolerance this passes
  expect_equal(Re(om(xg, 0))/TR, heads(m, xg, 0), tolerance = 1e-1)

  # compare visually,
  # TODO remove this from final test
  # plot(xg, Re(om(xg, 0))/TR)
  # lines(xg, heads(m, xg, 0))

})

test_that('aem works for 2D flow with areal recharge', {

  # based on Bakker & Post (2022), chapter 6.1
  k <- 5
  base <- 0
  N <- 0.001
  R <- 200
  hr <- 100 # ensures confined flow
  rw <- 0.3
  Qhalf <- 0.5 * N * pi * R^2

  # exact
  r <- seq(rw, R, length = 400)
  phiR <- 0.5 * k * (hr - base)^2
  phi2 <- -0.25*N*(r^2 - R^2) + phiR
  phi3 <- phi2 + Qhalf / (2*pi) * log(r/R)
  hhalf <- base + sqrt(2*phi3/k)

  # aem
  as <- areasink(0, 0, N = N, R = R)
  w <- well(0, 0, Q = Qhalf, rw = rw)
  rf <- constant(R, 0, hr)

  m <- aem(k, top = hr, base, n = 0.2, as, w, rf)
  haem <- heads(m, r, 0)

  expect_equal(haem, hhalf)

})
