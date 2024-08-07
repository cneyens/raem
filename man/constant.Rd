% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constant.R
\name{constant}
\alias{constant}
\title{Create a constant-head analytic element}
\usage{
constant(xc, yc, hc, ...)
}
\arguments{
\item{xc}{numeric, x location of the reference point.}

\item{yc}{numeric, y location of the reference point.}

\item{hc}{numeric, hydraulic head at the reference point.}

\item{...}{ignored}
}
\value{
Constant-head analytic element point which is an object of class \code{constant} and inherits from \code{element}.
}
\description{
\code{\link[=constant]{constant()}} creates an analytic element containing a constant head, often referred to
as \emph{the reference point}.
}
\examples{
rf <- constant(xc = -100, yc = 0, hc = 10)

}
