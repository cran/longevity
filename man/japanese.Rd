% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{japanese}
\alias{japanese}
\title{Japanese survival data}
\format{
A data frame with 1038 rows and 4 variables:
\describe{
\item{age}{integer, age (to the smallest year) at death (in years)}
\item{byear}{integer, birth year}
\item{count}{integer, number of death for cohort at given age}
\item{gender}{factor, the gender of the individuals; either \code{male} or \code{female}}
}
}
\source{
Table extracted from Hanayama & Sibuya (2016).
}
\usage{
japanese
}
\description{
This data frame contains information about the counts
of dead Japanese by gender and year of birth (cohort), categorized
by the whole part of age attained at death.
}
\details{
These data were obtained from the Annual Vital Statistics Report of Japan, released by the
Japanese government every year since 1947. The authors note that
"All the members of that cohort have died by the end of the observation period,
a procedure referred to as the extinct cohort method".
The data were obtained from the Human Mortality Database by the authors.
Only positive counts are reported and two records (Misao Okawa and Jiroemon Kimura) are
excluded because they do not correspond to the same selection mechanism.
}
\references{
Hanayama, N. and M. Sibuya (2016). Estimating the Upper Limit of Lifetime Probability Distribution, Based on Data of Japanese Centenarians, \emph{The Journals of Gerontology: Series A}, 71(\bold{8}), 1014–1021. \doi{10.1093/gerona/glv113}
}
\keyword{datasets}
