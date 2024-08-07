% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aem.R
\name{add_element}
\alias{add_element}
\alias{remove_element}
\title{Add or remove an element from an existing \code{aem} object}
\usage{
add_element(aem, element, name = NULL, solve = FALSE, ...)

remove_element(aem, name = NULL, type = NULL, solve = FALSE, ...)
}
\arguments{
\item{aem}{\code{aem} object.}

\item{element}{analytic element of class \code{element}.}

\item{name}{optional name of the element as character. Duplicate element names in \code{aem} are not allowed..}

\item{solve}{logical, should the model be solved after adding or removing the element? Defaults to \code{FALSE}.}

\item{...}{ignored}

\item{type}{class of the element(s) to remove. Either \code{name} or \code{type} should be specified in \code{\link[=remove_element]{remove_element()}}.}
}
\value{
The \code{aem} model with the addition of \code{element} or with the removal of element(s). If \code{solve = TRUE},
the model is solved using \code{\link[=solve.aem]{solve.aem()}}. The name of the new element is taken from the \code{name} argument,
the object name or set to \code{element_1} with \code{1} being the index of the new element in the element list. See examples.
}
\description{
\code{\link[=add_element]{add_element()}} adds a new element to an \code{aem} object.

\code{\link[=remove_element]{remove_element()}} removes an element from an \code{aem} object based on its name or type.
}
\examples{
m <- aem(k = 10, top = 10, base = 0, n = 0.2)
mnew <- add_element(m, constant(xc = 0, yc = 1000, hc = 12), name = 'rf')

# if name not supplied, tries to obtain it from object name
rf <- constant(xc = 0, yc = 1000, hc = 12)
mnew <- add_element(m, rf)

# or else sets it sequentially from number of elements
mnew <- add_element(m, constant(xc = 0, yc = 1000, hc = 12))

\dontshow{if (getRversion() >= '4.1.0') (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# add_element() adn remove_element() are pipe-friendly
mnew <- aem(k = 10, top = 10, base = 0, n = 0.2) |>
    add_element(rf, name = 'rf') |>
    add_element(headwell(xw = 0, yw = 100, rw = 0.3, hc = 8),
                name = 'headwell', solve = TRUE)
\dontshow{\}) # examplesIf}

# removing elements
mnew <- remove_element(mnew, name = 'rf')
mnew <- remove_element(mnew, type = 'headwell')

}
\seealso{
\code{\link[=aem]{aem()}}
}
