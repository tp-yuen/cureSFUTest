% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{grenander.test}
\alias{grenander.test}
\title{Testing sufficient follow-up using Grenander Estimator \eqn{\hat{f}^{G}}}
\usage{
grenander.test(y, delta, tau.c = NULL)
}
\arguments{
\item{y}{A numeric vector of the observed survival times, \eqn{Y_i}.}

\item{delta}{A numeric vector of the censoring indicators, \eqn{\Delta_i}.}

\item{tau.c}{A numeric of the end of the study time.
\code{NULL} indicates the maximum observed survival time \code{max(y)} is used.}
}
\value{
A list containing the results including the estimate
\eqn{\hat{f}^{G}(\tau_{G})}.
}
\description{
Testing the hypotheses:
\deqn{
   \tilde{H}_{0}: q_{1 - \epsilon} \geq \tau_{G}
   \quad\text{versus}\quad
   \tilde{H}_{a}: q_{1 - \epsilon} < \tau_{G},
}
using the Grenander estimator, \eqn{\hat{f}^{G}}.
}
\keyword{internal}
