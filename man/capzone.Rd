% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracelines.R
\name{capzone}
\alias{capzone}
\title{Calculate capture zone of a well element}
\usage{
capzone(aem, well, time, npar = 30, dt = time/20, zstart = aem$base, ...)
}
\arguments{
\item{aem}{\code{aem} object}

\item{well}{well analytic element of class \code{well} or inherits from it.}

\item{time}{numeric, time of the capture zone}

\item{npar}{integer, number of particles to use in the backward tracking. Defaults to 30.}

\item{dt}{numeric, time step length used in the particle tracking. Defaults \code{time / 20}.}

\item{zstart}{numeric value with the starting elevation of the particles. Defaults to the base of the aquifer.}

\item{...}{additional arguments passed to \code{\link[=tracelines]{tracelines()}}.}
}
\value{
\code{\link[=capzone]{capzone()}} returns an object of class \code{tracelines}.
}
\description{
\code{\link[=capzone]{capzone()}} determines the capture zone of a well element in the flow field by performing backward
particle tracking until the requested time is reached.
}
\details{
\code{\link[=capzone]{capzone()}} is a thin wrapper around \code{\link[=tracelines]{tracelines()}}. Backward particle tracking is performed using \code{\link[=tracelines]{tracelines()}}
and setting \code{forward = FALSE}. Initial particle locations are computed by equally spacing \code{npar} locations at the well
radius at the \code{zstart} elevation. To obtain a sharper delineation of the capture zone, try using more particles or
decreasing \code{dt}.

Note that different \code{zstart} values only have an effect in models with vertical flow components.
}
\examples{

k <- 10
top <- 10; base <- 0
n <- 0.3

uf <- uniformflow(TR = 100, gradient = 0.001, angle = -10)
rf <- constant(TR, xc = -1000, yc = 0, hc = 20)
w1 = well(200, 50, Q = 250)
w2 = well(-200, -100, Q = 450)

m <- aem(k, top, base, n = n, uf, rf, w1, w2)

cp5 <- capzone(m, w1, time = 5*365)
cp10 <- capzone(m, w2, time = 10*365)

xg <- seq(-800, 800, length = 100)
yg <- seq(-500, 500, length = 100)

contours(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
plot(cp5, add = TRUE)
plot(cp10, add = TRUE, col = 'orange3')

# model with vertical flow components
as <- areasink(0, 0, N = 0.001, R = 1500)
m <- aem(k, top, base, n = n, uf, rf, w1, w2, as)

# two different starting levels
cp5a <- capzone(m, w1, time = 5*365, zstart = base)
cp5b <- capzone(m, w1, time = 5*365, zstart = 8)

contours(m, xg, yg, col = 'dodgerblue3', nlevels = 20)
plot(cp5a, add = TRUE)
plot(cp5b, add = TRUE, col = 'forestgreen') # smaller zone

# plot the convex hull of the endpoints as a polygon
endp <- endpoints(cp5b)
hull <- chull(endp[, c('x', 'y')])
polygon(endp[hull, c('x', 'y')], col = adjustcolor('forestgreen', alpha.f = 0.7))

}
\seealso{
\code{\link[=tracelines]{tracelines()}}
}
