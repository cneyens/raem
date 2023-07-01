
#
# # examples
# TR = 100
# w <- well(xw = 50, yw = 0, Q = 200)
# wi <- well(xw = 0, yw = 0, Q = -100)
# uf <- uniformflow(gradient = 0.002, angle = -45, TR = TR)
# ml <- aem(TR = 100, w, wi, uf)
#
# xg <- seq(-100, 100, length = 500)
# yg <- seq(-75, 75, length = 100)
#
# contour(ml, xg, yg, 'potential', nlevels = 20, col = 'dodgerblue3')
# contour(ml, xg, yg, 'stream', nlevels = 20, col = 'orange3', drawlabels = FALSE, add = TRUE)
# plot(ml, add = TRUE)
#
# ml <- aem(TR = 100)
# w1 <- well(xw = 0, yw = 0, Q = 100, rw = 0.3)
# w2 <- headwell(xw = 400, yw = 300, rw = 0.3, hw = 20, TR = 100)
# uf <- uniformflow(TR = 100, gradient = 0.002, angle = 0)
# rf <- constant(TR = 100, xc = 400, yc = 0, hc = 22)
# ml <- aem(TR = 100, w1, w2, uf, rf)
#
# xg <- seq(-200, 600, length = 100)
# yg <- seq(-100, 500, length = 100)
# contour(ml, xg, yg, levels = seq(20, 40, 0.2), col = 'dodgerblue3', xlab = 'x (m)', ylab = 'y (m)', grid = grid())
# contour(ml, xg, yg, 'stream', levels = seq(-200, 100, ml$elements$w2$parameter/20), col = 'orange3', add = TRUE, drawlabels = FALSE)
# plot(ml, add = TRUE)
#
# ls1 <- linesink(x0 = -200, y0 = -150, x1 = 200, y1 = 150, sigma = 0.1)
# ml <- aem(TR = 100, ls1)
#
# xg <- seq(-400, 400, length = 100)
# yg <- seq(-300, 300, length = 100)
# h <- head(ml, xg, yg, as.grid = TRUE)
# xs <- seq(-200, 200, length = 100)
# ys <- seq(-150, 150, length = 100)
# hs <- head(ml, xs, ys)
#
# contour(ml, xg, yg, col = 'dodgerblue3', xlab = 'x (m)', ylab = 'y (m)', grid = grid())
# contour(ml, xg, yg, 'stream', col = 'orange3', add = TRUE, drawlabels = FALSE)
# plot(ml, add = TRUE)
# plot(sqrt((xs+200)^2 + (ys+150)^2), hs, type = 'l')
#
# TR <- 100
# hriver <- 10
# ml <- aem(TR = TR)
# ml <- add_element(ml, constant(TR = TR, xc = 0, yc = 1000, hc = hriver + 2), 'rf')
# ml <- add_element(ml, headwell(TR = TR, xw = 0, yw = 100, rw = 0.3, hw = hriver - 2), 'headwell', solve = FALSE)
# xls <- seq(-1600, 1600, length = 101)
# xls <- c(seq(-1600, -400, 200), seq(-350, 400, 50), seq(600, 1601, 200))
# yls <- 50*sin(pi*xls/400)
# for(i in seq_len(length(xls)-1)) {
#   hls <- headlinesink(TR = TR, x0 = xls[[i]], y0 = yls[[i]], x1 = xls[[i+1]], y1 = yls[[i + 1]], hc = hriver)
#   ml <- add_element(ml, hls, paste("headlinesink", i, sep = '_'), solve = FALSE)
# }
#
# ml <- solve(ml)
#
# xg1 <- seq(-1800, 1800, length = 100)
# yg1 <- seq(-1200, 1200, length = 101)
# xg2 <- seq(-400, 400, length = 100)
# yg2 <- seq(-100, 400, length = 100)
# h1 <- head(ml, xg1, yg1, as.grid = TRUE)
# h2 <- head(ml, xg2, yg2, as.grid = TRUE)
#
# contour(ml, xg1, yg1, nlevels = 10, col = 'dodgerblue3')
# plot(ml, add = TRUE)
# contour(ml, xg2, yg2, nlevels = 20, col = 'dodgerblue3')
# plot(ml, add = TRUE)
# contour(ml, xg2, yg2, z = 'stream', col = 'orange3', add = TRUE)
#
# TR = 100
# rf <- constant(TR = TR, xc = 0, yc = 0, hc = 20)
# as1 <- areasink(xc = -500, yc = 0, N = 0.001, R=500)
# as2 <- areasink(xc = 500, yc = 0, N = -0.001, R = 500)
# ml <- aem(TR = TR, rf, as1, as2)
# xg <- seq(-1500, 1500, length = 100)
# yg <- seq(-800, 800, length = 101)
# h <- head(ml, xg, yg, as.grid = T)
#
# contour(ml, xg, yg, nlevels = 10, col = 'dodgerblue3')
# filled.contour(xg, yg, t(h[dim(h)[1]:1,]), color.palette = function(n) hcl.colors(n, 'viridis'))
# plot(xg, h[,51], type = 'l')
#
# Q <- disvec(ml, xg, y = 100)
# plot(xg, Q[,1], type = 'l')
# lines(xg, Q[,2], type = 'l', col = 'red')
#
# hds <- head(ml, c(-450, -450.00001), y = 0)
# grad <- diff(hds)/0.00001
# Qx <- disvec(ml, -450, y = 0)[,1]
#
# xls0 <- c(0, 100, 200, 400, 600, 800, 1000, 1100, 1200)
# yls0 <- c(200, 200, 100, 100, 0, 0, 100, 300, 450)
# hls0 <- seq(39, 40.4, length = 8)
# xls1 <- c(0, 0, 200, 400, 600, 800, 1000, 1100, 1200)
# yls1 <- c(200, 400, 600, 600, 700, 700, 750, 800, 850)
# hls1 <- seq(39, 40.4, length = 8)
#
# TR <- 100
# ml <- aem(TR) |>
#   add_element(constant(TR, xc = 0, yc = 800, hc = 39.5), 'rf') |>
#   add_element(well(xw = 500, yw = 250, Q = 100), 'w0') |>
#   add_element(well(xw = 800, yw = 500, Q = 100), 'w1') |>
#   add_element(areasink(xc = 600, yc= 400, N = 0.001, R = 700))
# for(i in seq_along(hls0)) {
#   ml <- add_element(ml, headlinesink(TR, x0 = xls0[i], y0 = yls0[i], x1 = xls0[i+1], y1 = yls0[i+1], hc = hls0[i]),
#                     name = paste('hls0', i, sep = '_'))
# }
# for(i in seq_along(hls1)) {
#   ml <- add_element(ml, headlinesink(TR, x0 = xls1[i], y0 = yls1[i], x1 = xls1[i+1], y1 = yls1[i+1], hc = hls1[i]),
#                     name = paste('hls1', i, sep = '_'))
# }
# ml <- solve(ml)
#
# xg <- seq(-100, 1300, length = 100)
# yg <- seq(-100, 900, length = 100)
# h <- head(ml, xg, yg, as.grid = TRUE)
#
# contour(ml, xg, yg, levels = seq(38, 41, 0.1), asp = 1, col = 'dodgerblue3', panel.first = grid())
# plot(ml, add = TRUE)
#
# filled.contour(xg, yg, h, color.palette = function(n) hcl.colors(n, 'viridis'),
#                plot.axes = {plot(ml, add = TRUE)})
#
# plot(ml, xlim = range(xg), ylim = range(yg), panel.first = grid(), asp = 1)
#
# h <- head(ml, xg, y = 400)
# plot(xg, h, type = 'l')
#
# TR = 100
# b = 10 # thickness
# w <- well(xw = 100, yw = 150, Q = 200)
# uf <- uniformflow(gradient = 0.002, angle = 50, TR = TR)
# rf <- constant(TR, -100, -1000, 10)
# ml <- aem(TR = 100, uf, rf, w)
#
# xg <- seq(0, 300, length = 100)
# yg <- seq(0, 300, length = 100)
# contour(ml, xg, yg)
#
# Q <- disvec(ml, xg, yg, as.grid = TRUE, magnitude = TRUE)
# contour(xg, yg, Q[,,3], nlevels = 100)
