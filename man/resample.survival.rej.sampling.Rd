% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{resample.survival.rej.sampling}
\alias{resample.survival.rej.sampling}
\title{Generate bootstrap sample of event time
using acceptance-rejection sampling}
\usage{
resample.survival.rej.sampling(nr, gcmlcm.fit, tau.c, n, skm.bw.multi = 1)
}
\arguments{
\item{nr}{An integer of number of bootstrap samples}

\item{gcmlcm.fit}{A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.}

\item{tau.c}{A numeric of the end of the study time.}

\item{n}{An integer of sample size.}

\item{skm.bw.multi}{A numeric of scaling to the bandwidth for smoothing
the density estimator for generating bootstrap samples.}
}
\value{
A numeric vector of bootstrap sample of event time.
}
\description{
Generate bootstrap sample of event time
using acceptance-rejection sampling
}
\keyword{internal}
