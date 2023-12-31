% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.R
\name{ks_test}
\alias{ks_test}
\title{Goodness-of-fit diagnostics}
\usage{
ks_test(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  ltrunc = NULL,
  rtrunc = NULL,
  type = c("right", "left", "interval", "interval2"),
  family = c("exp", "gp", "gomp", "gompmake", "weibull", "extgp", "gppiece",
    "extweibull", "perks", "beard", "perksmake", "beardmake"),
  B = 999L,
  arguments = NULL,
  ...
)
}
\arguments{
\item{time}{excess time of the event of follow-up time, depending on the value of event}

\item{time2}{ending excess time of the interval for interval censored data only.}

\item{event}{status indicator, normally 0=alive, 1=dead. Other choices are \code{TRUE}/\code{FALSE} (\code{TRUE} for death).
For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have experienced an event.}

\item{thresh}{vector of thresholds}

\item{ltrunc}{lower truncation limit, default to \code{NULL}}

\item{rtrunc}{upper truncation limit, default to \code{NULL}}

\item{type}{character string specifying the type of censoring. Possible values are "\code{right}", "\code{left}", "\code{interval}", "\code{interval2}".}

\item{family}{string; choice of parametric family}

\item{B}{number of bootstrap simulations}

\item{arguments}{a named list specifying default arguments of the function that are common to all \code{elife} calls}

\item{...}{additional parameters, currently ignored}
}
\value{
a list with elements
\itemize{
\item \code{stat} the value of the test statistic
\item \code{pval} p-value obtained via simulation
}
}
\description{
Warning: EXPERIMENTAL
Compute the Kolmogorov-Smirnov
test statistic and compare it with a simulated null
distribution obtained via a parametric bootstrap.
}
\note{
The bootstrap scheme requires simulating new data,
fitting a parametric model and estimating the nonparametric
maximum likelihood estimate for each new sample.
This is computationally intensive in large samples.
}
\keyword{internal}
