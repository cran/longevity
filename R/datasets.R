#' Dutch survival data
#'
#' This data frame contains information about all Dutch who died
#' above age 92 years between 1986 and 2015. Observations are
#' doubly truncated and such bounds are calculated based on the
#' range of plausible values for these variables.
#' There are 226 records that are interval-censored and interval-truncated
#' for which \code{bdate}, \code{ddate} and \code{ndays} is missing (\code{NA}).
#'
#' @format A data frame with 305143 rows and 11 variables:
#' \describe{
#' \item{ndays}{survival time (in days)}
#' \item{bdate}{the smallest plausible birth date given information
#' about month of birth and death and survival (\code{Date})}
#' \item{bmonth}{month of birth}
#' \item{byear}{year of birth}
#' \item{ddate}{the largest plausible death date given information
#'  about month of birth and death and survival (\code{Date})}
#' \item{dmonth}{month of death}
#' \item{dyear}{year of death}
#' \item{ltrunc}{minimum age (in days); the maximum of either 92
#'  years or the number of days reached in 1986}
#' \item{rtrunc}{maximum age (in days) an individual could have
#' reached by the end of 2015}
#' \item{gender}{factor indicating gender of individual, either
#' \code{female} or \code{male}}
#' \item{valid}{quality flag; \code{A} for individuals born in
#' the Netherlands, \code{B} for individuals born abroad who died
#'  in the Netherlands}
#' }
#' @references Einmahl, J.J., J.H.J. Einmahl and L. de Haan (2019). \emph{Limits to Human Life Span Through Extreme Value Theory}, Journal of the American Statistical Association, \bold{114}(527), 1075-1080. \doi{10.1080/01621459.2018.1537912}
#'
#' @source Statistics Netherlands (CBS). Accessed via the Supplemental material of Einmahl, Einmahl and de Haan (2019)
"dutch"

#' England and Wales semi-supercentenarian
#'
#' This data frame contains information about 3866
#' Welsh and English who died at age ranging from 105 to 110
#' between 2000 and 2014 (except for two women
#' who died late in December 1999) and a subset of UK supercentenarians from
#' the IDL database (5 male, 80 female) who died during the same period.
#' All records for people who died at age 109 and all men,
#' plus a stratified sample of the women were validated
#' by the General Register Office (GRO). Observations are
#' doubly truncated.
#'
#' In the original data forwarded by the IDL staff,
#' there were 7 dubious records (missing birth day or month)
#' that were excluded. The referenced technical reports describes the validation
#' procedure in more details and includes (approximate) sampling weights for the validation
#' sample of women who died age 105-108.
#' @references Office for National Statistics (2016). Accuracy of official high-age
#' population estimates, in England and Wales: an evaluation. Technical report,
#' \url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/methodologies/accuracyofofficialhighagepopulationestimatesinenglandandwalesanevaluation}
#' @source Ngaire Coombs, Office for National Statistics (ONS)
#' @format A data frame with 3951 rows and 7 variables:
#' \describe{
#' \item{ndays}{survival time (in days)}
#' \item{bdate}{birth date (\code{Date})}
#' \item{ddate}{death date (\code{Date})}
#' \item{ltrunc}{minimum age (in days); the maximum of 38350 days (approximately 105 years)
#' or the number of days reached in 2000}
#' \item{rtrunc}{maximum age (in days) an individual could have reached by the end of 2014}
#' \item{gender}{factor indicating gender of individual, either \code{female} or \code{male}}
#' \item{valid}{quality flag; \code{A} for validated records, \code{B} for unchecked records}
#' }
"englandwales"


#' England and Wales simulated supercentenarian data
#'
#' This data frame contains information about 179 fake records mimicking
#' Welsh and English who died age 110 and above
#'
#' @format A data frame with 179 rows and 3 variables:
#' \describe{
#' \item{time}{survival time above 110 (in years)}
#' \item{ltrunc}{minimum age above 110 (in years), or zero};
#' \item{rtrunc}{maximum age (in years) an individual could have reached by the end of the time frame}
#' }
"ewsim"

#' Japanese survival data
#'
#' This data frame contains information about the counts
#' of dead Japanese by gender and year of birth (cohort), categorized
#' by the whole part of age attained at death.
#'
#' These data were obtained from the Annual Vital Statistics Report of Japan, released by the
#' Japanese government every year since 1947. The authors note that
#' "All the members of that cohort have died by the end of the observation period,
#' a procedure referred to as the extinct cohort method".
#' The data were obtained from the Human Mortality Database by the authors.
#' Only positive counts are reported and two records (Misao Okawa and Jiroemon Kimura) are
#' excluded because they do not correspond to the same selection mechanism.
#'
#' @format A data frame with 1038 rows and 4 variables:
#' \describe{
#' \item{age}{integer, age (to the smallest year) at death (in years)}
#' \item{byear}{integer, birth year}
#' \item{count}{integer, number of death for cohort at given age}
#' \item{gender}{factor, the gender of the individuals; either \code{male} or \code{female}}
#' }
#' @references Hanayama, N. and M. Sibuya (2016). Estimating the Upper Limit of Lifetime Probability Distribution, Based on Data of Japanese Centenarians, \emph{The Journals of Gerontology: Series A}, 71(\bold{8}), 1014–1021. \doi{10.1093/gerona/glv113}
#'
#' @source Table extracted from Hanayama & Sibuya (2016).
"japanese"


#' Japanese survival data (2)
#'
#' This data frame is extracted from Table 10.3 from Chapter 10, "Centenarians and Supercentenarians in Japan", in the Monograph Exceptional lifespans. The data were constructed by the extinct cohort method and are stratified by age cohort (five year group, except 1899-1900) and by sex. Note that the family registry system (KOSEKI), introduced in 1872, was standardized in 1886.
#'
#' @format A data frame with 216 rows and 4 variables:
#' \describe{
#' \item{age}{integer, age (to the smallest year) at death (in years)}
#' \item{bcohort}{factor, birth cohort}
#' \item{count}{integer, number of death for cohort at given age}
#' \item{gender}{factor, the gender of the individuals; either \code{male} or \code{female}}
#' }
#' @references
#' Saito, Yasuhiko and Futoshi Ishii, and Jean-Marie Robine (2021). \emph{Centenarians and Supercentenarians in Japan}. In \emph{Exceptional lifespans}, Maier, H., Jeune, B., Vaupel, J. W. (Eds.), Demographic research monographs 17 VII, pp. 125-145. Cham, Springer.
#'
#' @source Table 10.3
"japanese2"



#' Italian semi-supercentenarian
#'
#' This data frame contains information about 3836 Italians
#' individually validated survival lifetimes times in days
#' of all persons in Italy who were at least 105 years old
#' at some point in the period from 1 January 2009
#' to 31 December 2015.
#' Observations are left-truncated and right-censored.
#' These data are not publicly available, but can be purchased
#' from the Italian National Institute of Statistics by
#' registering at the Contact Center and mentioning the
#' Semi-supercentenarian Survey and Marco Marsili
#' as contact person.
#'
#' @references Istituto Nazionale di Statistica
#' @format A data frame with 3836 rows and 6 variables:
#' \describe{
#' \item{ndays}{survival time (in days)}
#' \item{bdate}{birth date (\code{Date})}
#' \item{ddate}{death date (\code{Date}), or \code{NA_Date_}
#' if the person is alive at the end of the sampling}
#' \item{ltrunc}{minimum age (in days); the maximum of 38351 days
#' (approximately 105 years) or the number of days reached in 2009}
#' \item{event}{integer indicating the censoring pattern; \code{0} for
#' right-censored records, \code{1} for fully observed}
#' \item{gender}{factor indicating gender of individual,
#' either \code{female} or \code{male}}
#' }
"italian"


#' French semi-supercentenarian
#'
#' This data frame contains information about 9853 French semi-supercentenarian,
#' part of the International Database on Longevity (IDL). All
#' supercentenarian records were validated, but only a random sample
#' of semi-supercentenarians were validated.
#' Lifetimes are interval truncated;only people above 110 born
#' after 1978 and people above 105 born after 1987 are included.
#'
#' @references International Database on Longevity
#' @format A data frame with 9853 rows and 6 variables:
#' \describe{
#' \item{ndays}{survival time (in days)}
#' \item{bdate}{birth date (\code{Date})}
#' \item{ddate}{death date (\code{Date}) if the person
#' is alive at the end of the sampling}
#' \item{ltrunc}{minimum age (in days); the maximum of 38350
#'  days (approximately 105 years) or the number of days
#'  reached on January 1st 1978 (supercentenarian) or 1987
#'  (semisupercentenarian)}
#' \item{rtrunc}{maximum age (in days) an individual could
#' have reached by the end of 2017}
#' \item{gender}{factor indicating gender of individual,
#' either \code{female} or \code{male}}
#' }
"french"


#' International Database on Longevity (2021)
#'
#' This database contains data downloaded from \url{supercentenarian.org},
#' including third first waves, new data provided by ONS for semisupercentenarian,
#' data for Switzerland and Italy previously available for download and removed
#' for confidentiality. Data from Japan and from people aged less than 110 from the USA are excluded
#' because they are of dubious quality. For the USA, the semisupercentenarian records are validated,
#' but this is only a fraction of a cohort whose size is unknown and they are not representative of the whole population.
#' The birth and death dates of the USA people are unknown (only years are given,
#' so the largest plausible range is recorded given the survival in years).
#'
#' Only dead individuals are included, so the records are truncated.
#' For countries with semisupercentenarians and with different collection period for
#' semisupercentenarians (105-109) and supercentenarians (110+), there
#' are some configurations leading to double interval truncation, in which case
#' data are defined in \eqn{[\code{ltrunc1}, \code{rtrunc1}] \cup [\code{ltrunc2},\code{rtrunc2}]}.
#'
#' @references International Database on Longevity
#' @format A data frame with 17721 rows and 10 variables:
#' \describe{
#' \item{country}{factor, one of \code{CH} (Switzerland), \code{OS} (Austria), \code{BE} (Belgium), \code{QC} (Quebec), \code{DE} (Germany), \code{DN} (Denmark), \code{ES} (Spain), \code{FI} (Finland), \code{FR} (France), \code{NO} (Norway), \code{SV} (Sweden), \code{EW} (England and Wales), \code{IT} (Italy) and \code{US} (United States of America)}
#' \item{ndays}{ integer; survival (in days)}
#' \item{ageyear}{ integer; floor of maximum age (in years) reached at death}
#' \item{gender}{factor; \code{male} or \code{female}}
#' \item{bdate}{Date; birth date (except for US)}
#' \item{ddate}{Date; death date (except for US)}
#' \item{ltrunc1}{integer; lower truncation limit (in days); the minimum number of days someone would have survived to be included in the sampling frame (first interval)}
#' \item{rtrunc1}{integer; upper truncation limit (in days); the maximum number of days someone would have survived to be included in the sampling frame (first interval)}
#' \item{ltrunc2}{integer; lower truncation limit (in days); the minimum number of days someone would have survived to be included in the sampling frame (second interval) if applicable, \code{NA} otherwise}
#' \item{rtrunc2}{integer; upper truncation limit (in days); the maximum number of days someone would have survived to be included in the sampling frame (second interval) if applicable, \code{NA} otherwise}
#' }
"idl"


#' IDL metadata
#'
#' This data frame contains country codes
#' and the associated data collection period
#' corresponding to the range for age at death.
#'
#' @details Due to confidentiality restrictions, some data that were available in previous versions of the IDL for Switzerland, Italy and some entries for Japan and Belgium have been removed.
#' As the IDL metadata are updated somewhat regularly and former versions of the database are not preserved, results from published analyses are replicable but not reproducible.
#'
#' @references International Database on Longevity, extracted on February 13th, 2023
#' @format A data frame with 21 rows and 4 variables:
#' \describe{
#' \item{country}{factor, one of \code{AUT} (Austria), \code{BEL} (Belgium), \code{CAN} (Quebec), \code{DEU} (Germany), \code{DNK} (Denmark), \code{ESP} (Spain), \code{FIN} (Finland), \code{FRA} (France), \code{JPN} (Japan), \code{NOR} (Norway), \code{SWE} (Sweden), \code{EAW} (England and Wales) and \code{USA} (United States of America)}
#' \item{group}{factor, either \code{105-109} for semi-supercentenarians or \code{110+} for supercentenarians"}
#' \item{ldate}{Date, smallest death date}
#' \item{rdate}{Date, latest death date}
#' }
"idlmetadata"


#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics hist
#' @importFrom graphics mtext
#' @importFrom graphics rug
#' @importFrom stats approxfun
#' @importFrom stats dbinom
#' @importFrom stats pbinom
#' @importFrom stats qbinom
#' @importFrom stats dexp
#' @importFrom stats dweibull
#' @importFrom stats family
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats pchisq
#' @importFrom stats pexp
#' @importFrom stats ppoints
#' @importFrom stats predict
#' @importFrom stats pweibull
#' @importFrom stats qchisq
#' @importFrom stats qexp
#' @importFrom stats qnorm
#' @importFrom stats qweibull
#' @importFrom stats rexp
#' @importFrom stats rmultinom
#' @importFrom stats runif
#' @importFrom stats rweibull
#' @importFrom stats smooth.spline
#' @importFrom stats weighted.mean
NULL
