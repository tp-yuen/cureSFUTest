% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothing.R
\name{boundary.smoother}
\alias{boundary.smoother}
\title{Kernel smoothing with boundary correction using tri-weight kernel}
\usage{
boundary.smoother(x, jumptimes, jumpheights, bandwidth, a, b)
}
\arguments{
\item{x}{A numeric vector of points at which the smoothing is computed.}

\item{jumptimes}{A numeric vector of points at which there is a jump.}

\item{jumpheights}{A numeric vector of jump heights.}

\item{bandwidth}{A positive number of the bandwidth parameter.}

\item{a}{A numeric of left endpoint of the support.}

\item{b}{A numeric of right endpoint of the support.}
}
\value{
A numeric vector of the computed values at \code{x}.
}
\description{
Kernel smoothing with boundary correction using tri-weight kernel
}
\keyword{internal}
