% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.censTruncLimits}
\alias{.censTruncLimits}
\title{Identification sets}
\usage{
.censTruncLimits(tsets, lcens, rcens, ltrunc, rtrunc, trunc, cens)
}
\arguments{
\item{tsets}{Turnbull's sets}

\item{lcens}{numeric vector of left censoring}

\item{rcens}{numeric vector of right censoring}

\item{ltrunc}{numeric vector of left truncation}

\item{rtrunc}{numeric vector of right truncation}

\item{trunc}{logical are observation truncated?}
}
\value{
a matrix with the bounds of the intervals for Turnbull sets
}
\description{
Identification sets
}
\keyword{internal}
