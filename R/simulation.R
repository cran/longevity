#' Simulate excess lifetime with truncation or right-censoring
#'
#' This function dispatches simulations accounting for potential left-truncation (remove by setting lower to zero).
#' If \code{type2=ltrt}, simulated observations will be lower than the upper bounds \code{upper}.
#' If \code{type2=ltrc}, simulated observations are capped at \code{upper} and the observation is right-censored (\code{rcens=TRUE}).
#'
#' @param n sample size
#' @param scale scale parameter(s)
#' @param rate rate parameter(s)
#' @param shape shape parameter(s)
#' @param family string; choice of parametric family
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds
#' @param type2 string, either \code{none}, \code{ltrt} for left- and right-truncated data or \code{ltrc} for left-truncated right-censored data
#' @export
#' @return either a vector of observations or, if \code{type2=ltrc}, a list with \code{n} observations \code{dat} and a logical vector of the same length with \code{TRUE} for right-censored observations and \code{FALSE} otherwise.
#' @examples
#' set.seed(1234)
#' n <- 500L
#' # Simulate interval truncated data
#' x <- samp_elife(n = n,
#'                 scale = 2,
#'                 shape = 1.5,
#'                 lower = low <- runif(n),
#'                 upper = upp <- runif(n, min = 3, max = 15),
#'                 type2 = "ltrt",
#'                 family = "weibull")
#' coef(fit_elife(
#'    time = x,
#'    ltrunc = low,
#'    rtrunc = upp,
#'    family = "weibull"))
#' # Simulate left-truncated right-censored data
#' x <- samp_elife(n = n,
#'                 scale = 2,
#'                 shape = 1.5,
#'                 lower = low <- runif(n),
#'                 upper = upp <- runif(n, min = 3, max = 15),
#'                 type2 = "ltrc",
#'                 family = "gomp")
#' #note that the return value is a list...
#' coef(fit_elife(
#'    time = x$dat,
#'    ltrunc = low,
#'    event = !x$rcens,
#'    family = "gomp"))
samp_elife <- function(n,
                       scale,
                       rate,
                       shape = NULL,
                       lower = 0,
                       upper = Inf,
                       family = c("exp",
                                  "gp",
                                  "gomp",
                                  "gompmake",
                                  "weibull",
                                  "extgp",
                                  "gppiece",
                                  "extweibull",
                                  "perks",
                                  "beard",
                                  "perksmake",
                                  "beardmake"),
                       type2 = c("none",
                                 "ltrt",
                                 "ltrc",
                                 "ditrunc")){
  family <- match.arg(family)
  type2 <- match.arg(type2)
  if(type2 == "none"){
    check_elife_dist(rate = rate, scale = scale, shape = shape, family = family)
    relife(n = n,
           rate = rate,
           scale = scale,
           shape = shape,
           family = family)
  } else if(type2 == "ltrt"){
    stopifnot(length(lower) %in% c(1L, n),
              length(upper) %in% c(1L, n),
              isTRUE(all(lower < upper)))
    r_dtrunc_elife(n = n,
                   scale = scale,
                   rate = rate,
                   shape = shape,
                   lower = lower,
                   upper = upper,
                   family = family)
  } else if(type2 == "ltrc"){
    stopifnot(length(lower) %in% c(1L, n),
              length(upper) %in% c(1L, n),
              isTRUE(all(lower < upper)))
    r_ltrc_elife(n = n,
                 scale = scale,
                 rate = rate,
                 shape = shape,
                 lower = lower,
                 upper = upper,
                 family = family)
  } else if(type2 == "ditrunc"){
    r_ditrunc_elife(n = n,
                    rate = rate,
                    scale = scale,
                    shape = shape,
                    lower = lower,
                    upper = upper,
                    family = family)
  }
}



#' Sample observations from an interval truncated excess lifetime distribution
#'
#' @inheritParams samp_elife
#' @return a vector of \code{n} observations
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_dtrunc_elife(n = n, scale = 1, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = runif(n, 6, 10),
#'                family = "gp")
r_dtrunc_elife <- function(n,
                          scale,
                          rate,
                          shape,
                          lower,
                          upper,
                          family = c("exp",
                                     "gp",
                                     "gomp",
                                     "gompmake",
                                     "weibull",
                                     "extgp",
                                     "gppiece",
                                     "extweibull",
                                     "perks",
                                     "beard",
                                     "perksmake",
                                     "beardmake")
                          ){
  family <- match.arg(family)
  check_elife_dist(rate = rate, scale = scale, shape = shape, family = family)
  if(length(lower) > 1){
    stopifnot("Sample size should match length of 'lower'" = length(lower) == n)
  }
  if(length(upper) > 1){
    stopifnot("Sample size should match length of 'upper'" = length(upper) == n)
  }
  samp <- qelife(
    p = pelife(q = lower, scale = scale, shape = shape, rate = rate, family = family) +
      runif(n) * (
        pelife(q = upper, scale = scale, shape = shape, rate = rate, family = family) -
          pelife(q = lower, scale = scale, shape = shape, rate = rate, family = family)),
    scale = scale,
    shape = shape,
    rate = rate,
    family = family)
  if(isTRUE(any(!is.finite(samp)))){
    stop("Inversion method failed: infinite values or NAs generated")
  }
  return(samp)
}



#' Sample observations from an interval truncated excess lifetime distribution
#'
#' @inheritParams samp_elife
#' @return a vector of \code{n} observations
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_dtrunc_elife(n = n, scale = 1, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = runif(n, 6, 10),
#'                family = "gp")
r_ditrunc_elife <- function(n,
                            rate,
                           scale,
                           shape,
                           lower,
                           upper,
                           family = c("exp",
                                      "gp",
                                      "gomp",
                                      "gompmake",
                                      "weibull",
                                      "extgp",
                                      "gppiece",
                                      "extweibull",
                                      "perks",
                                      "beard",
                                      "perksmake",
                                      "beardmake")

){
  stopifnot("`lower` should be a matrix." = is.matrix(lower),
            "`upper` should be a matrix." = is.matrix(upper),
            "`lower` must be a matrix with 2 columns." = ncol(lower) == 2,
            "`upper` must be a matrix with 2 columns." = ncol(upper) == 2,
            "`lower` and `upper` should have the same number of rows." = nrow(lower) == nrow(upper),
            "`lower` and `upper` should have 1 or n rows." = nrow(upper) == 1 | nrow(upper) == n)
  family <- match.arg(family)
  check_elife_dist(rate = rate, scale = scale, shape = shape, family = family)
  probint1 <- pelife(q = upper[,1], rate = rate, scale = scale, shape = shape, family = family) -
    pelife(q = lower[,1], rate = rate, scale = scale, shape = shape, family = family)
  probint2 <- ifelse(is.na(lower[,2]), 0,
                     pelife(q = upper[,2], rate = rate, scale = scale, shape = shape, family = family) -
                       pelife(q = lower[,2], rate = rate, scale = scale, shape = shape, family = family)
  )
  if(isTRUE(any((probint1 + probint2) > 1))){
    stop("Invalid input")
  }

  cutoff <- probint1 / (probint1 + probint2)
  unif <- runif(n)
  qelife(ifelse(unif < cutoff,
         unif*(probint1 + probint2) + pelife(q = lower[,1], rate = rate, scale = scale, shape = shape, family = family),
         unif*(probint1 + probint2) + pelife(q = lower[,2], rate = rate, scale = scale, shape = shape, family = family) - probint1),
         rate = rate, scale = scale, shape = shape, family = family)
}


#' Sample observations from a left-truncated right-censored excess lifetime distribution
#'
#' @param n sample size
#' @param scale scale parameter
#' @param shape shape parameter(s)
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}) or extended generalized Pareto (\code{extgp}).
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds above which data are right-truncated
#' @return a list with \code{n} observations \code{dat} and a logical vector of the same length with \code{TRUE} for right-censored observations and \code{FALSE} otherwise.
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_ltrc_elife(n = n, scale = 5, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = 5,
#'                family = "gp")
r_ltrc_elife <- function(n,
                         scale,
                         rate,
                         shape,
                         lower,
                         upper,
                         family = c("exp",
                                    "gp",
                                    "gomp",
                                    "gompmake",
                                    "weibull",
                                    "extgp",
                                    "gppiece",
                                    "extweibull",
                                    "perks",
                                    "beard",
                                    "perksmake",
                                    "beardmake")

){
  family <- match.arg(family)
  check_elife_dist(rate = rate, scale = scale, shape = shape, family = family)
  if(length(lower) > 1){
    stopifnot("Sample size should match length of 'lower'" = length(lower) == n)
  }
  if(length(upper) > 1){
    stopifnot("Sample size should match length of 'upper'" = length(upper) == n)
  }
  samp <- qelife(p = pelife(lower, scale = scale, shape = shape, rate = rate, family = family) +
         runif(n)*pelife(lower, scale = scale, shape = shape, rate = rate, family = family, lower.tail = FALSE),
         scale = scale, shape = shape, rate = rate, family = family)
  rcens <- ifelse(samp < upper, FALSE, TRUE)
  samp <- pmin(samp, upper)
  if(any(is.infinite(samp))){
    warning("Infinite values returned by procedure due to numerical overflow.")
  }
  return(list(dat = samp,
              rcens = rcens)
         )
}

#' Sample lifetime from excess lifetime model
#'
#' Given parameters of a \code{elife} distribution, sampling window and
#' birth dates with excess lifetimes, sample new observations; excess lifetime
#' at \code{c1} are sampled from an exponential distribution, whereas
#' the birth dates are sampled from a jittered histogram-based distribution
#' The new excess lifetime above the threshold are right-censored if they exceed
#' \code{c2}.
#'
#' @inheritParams r_dtrunc_elife
#' @param xcal date at which individual reaches \code{u} years
#' @param c1 date, first day of the sampling frame
#' @param c2 date, last day of the sampling frame
#' @return list with new birthdates (\code{xcal}), excess lifetime at \code{c1} (\code{ltrunc}),
#' excess lifetime above \code{u} (\code{dat}) and right-censoring indicator (\code{rightcens}).
#' @export
#' @keywords internal
samp2_elife <- function(n,
                         scale,
                         shape,
                        family = c("exp",
                                   "gp",
                                   "gomp",
                                   "gompmake",
                                   "weibull",
                                   "extgp",
                                   "gppiece",
                                   "extweibull",
                                   "perks",
                                   "beard",
                                   "perksmake",
                                   "beardmake"),
                         xcal,
                         c1,
                         c2){
  xcal <- as.Date(xcal)
  stopifnot("Date `c1` is missing" = !missing(c1),
            "Date `c2` is missing" = !missing(c2),
            "`c1` should be a single date" = length(c1) == 1L,
            "`c2` should be a single date" = length(c2) == 1L)
  c1 <- as.Date(c1)
  c2 <- as.Date(c2)
  ltrunc <- as.numeric(pmax(0, c1 - xcal))
  sample_dates <- function(n, xcal, c1, c2, ltrunc){
    sample_ltrunc <- function(n, ltrunc){
      sort(round(rexp(n, rate = 1/mean(365.25*ltrunc[ltrunc>0])))/365.25, decreasing = TRUE)
    }
    nday <- as.numeric(xcal-c1)
    nday <- nday[nday>0]
    nmax <- as.numeric(c2-c1)
    ssltrunc <- sample_ltrunc(round(sum(ltrunc>0)*n/length(ltrunc)), ltrunc = ltrunc)
    xhist <- hist(nday, plot = FALSE)
    bins <- with(xhist, sample(length(mids), n-length(ssltrunc), p=density, replace=TRUE)) # choose a bin
    result <- round(runif(length(bins), xhist$breaks[bins], pmin(xhist$breaks[bins+1], nmax-1)))
    list(xcal = as.Date(round(c(-ssltrunc*365.25, sort(result))), origin = c1),
         ltrunc = c(ssltrunc, rep(0, n-length(ssltrunc))))
  }

  traject <- r_dtrunc_elife(n = n,
                           scale = scale,
                           shape = shape,
                           lower = 0,
                           upper = Inf,
                           family = family)
  sdates <- sample_dates(n = n, xcal, c1, c2, ltrunc)
  lifeah <- pmax(1,pmin(round(365.25*(traject + sdates$ltrunc)), as.numeric(c2-sdates$xcal)))
  rcens_new <- lifeah == as.numeric(c2-sdates$xcal)
  list(xcal = sdates$xcal, ltrunc = sdates$ltrunc, dat = lifeah/365.25, rightcens = rcens_new)
}
