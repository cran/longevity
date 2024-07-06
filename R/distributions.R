
#' Distribution function of the generalized Pareto distribution
#'
#' @param q vector of quantiles.
#' @param loc location parameter.
#' @param scale scale parameter, strictly positive.
#' @param shape shape parameter.
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname gpd
pgpd <- function(q,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE,
                 log.p = FALSE) {
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc))
  )
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
  q <- pmax(q - loc, 0) / scale
  if (shape == 0) {
    p <- 1 - exp(-q)
  } else {
    p <- 1 - exp((-1 / shape) * log1p(pmax(-1, shape * q)))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Density function of the generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname gpd
dgpd <- function (x,
                  loc = 0,
                  scale = 1,
                  shape = 0,
                  log = FALSE) {
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc))
  )
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    return(dexp(
      x = x - loc,
      rate = 1 / scale,
      log = log
    ))
  }
  d <- (x - loc) / scale
  index <- (d >= 0 & ((1 + shape * d) >= 0)) & !is.na(d)
  d[index] <- log(1 / scale) -
    (1 / shape + 1) * log(1 + shape * d[index])
  d[!index & !is.na(d)] <- -Inf
  d[is.na(d)] <- NA
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

#' Hazard function of the exponential distribution
#'
#' @param x vector of quantiles
#' @param rate rate vector of rates
#' @param log logical; if \code{FALSE} (default), return the log hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname exponential
hexp <- function(x, rate = 1, log = FALSE){
  stopifnot(isTRUE(all(rate > 0)))
  ifelse(x <= 0 | is.na(x), NA, rate)
}

#' Hazard function of the generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the log hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname gpd
hgpd <- function (x,
                  loc = 0,
                  scale = 1,
                  shape = 0,
                  log = FALSE) {
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc))
  )
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    return(hexp(
      x = x - loc,
      rate = 1 / scale
    ))
  }
  haz <-   ifelse(x >= 0 & ((shape * x/scale) >= -1) & !is.na(x),
                 -log(pmax(0, scale + x*shape)),
                 NA)
  if(log){
    return(haz)
  } else{
    return(exp(haz))
  }
}
#' Quantile function of the generalized Pareto distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @return a vector of quantiles
#' @export
#' @keywords internal
#' @rdname gpd
qgpd <- function(p,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("\"p' must contain probabilities in (0,1)")
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc))
  )
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
  if (lower.tail) {
    p <- 1 - p
  }
  if (shape == 0) {
    return(loc - scale * log(p))
  } else {
    return(loc + scale * (p ^ (-shape) - 1) / shape)
  }
}


#' Distribution function of the extended Weibull distribution
#'
#' @param q vector of quantiles.
#' @param scale scale parameter, strictly positive.
#' @param shape1 shape parameter of the generalized Pareto component.
#' @param shape2 shape parameter of the Weibull component.
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname extweibull
pextweibull <- function(q,
                        scale = 1,
                        shape1 = 0,
                        shape2 = 1,
                        lower.tail = TRUE,
                        log.p = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extweibull")
  q <- (pmax(0, q) / scale)^shape2
  if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))){
  # Weibull model
      p <- 1 - exp(-q)
  } else {
    p <- 1 - exp((-1 / shape1) * log1p(pmax(-1, shape1 * q)))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Density function of the extended Weibull distribution
#'
#' @inheritParams pextweibull
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname extweibull
dextweibull <- function (x,
                         scale = 1,
                         shape1 = 0,
                         shape2 = 1,
                         log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extweibull")
  if (isTRUE(all.equal(c(shape1, shape2), c(0, 1), check.attributes = FALSE))) {
    return(dexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  } else if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))) {
    return(dweibull(
      x = x,
      shape = shape2,
      scale = scale,
      log = log
    ))
  } else if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    return(dgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape1,
      log = log
    ))
  }
  d <- x / scale
  index <- (d >= 0 & ((1 + shape1 * d^shape2) >= 0)) & !is.na(d)
  d[index] <- -log(scale) + (shape2 - 1) * log(d[index]) - (1 / shape1 + 1) * log(1 + shape1 * d[index]^shape2) + log(shape2)
  d[!index & !is.na(d)] <- -Inf
  d[is.na(d)] <- NA
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

#' Hazard function of the extended Weibull distribution
#'
#' @inheritParams pextweibull
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the hazard, else the log hazard.
#' @return a vector of (log)-hazard
#' @export
#' @keywords internal
#' @rdname extweibull
hextweibull <- function (x,
                         scale = 1,
                         shape1 = 0,
                         shape2 = 1,
                         log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extweibull")
  if (isTRUE(all.equal(c(shape1, shape2), c(0, 1), check.attributes = FALSE))) {
    return(hexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  } else if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))) {
    return(hweibull(
      x = x,
      shape = shape2,
      scale = scale,
      log = log
    ))
  } else if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    return(hgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape1,
      log = log
    ))
  }
  haz <- suppressWarnings(ifelse(x >= 0 & ((1 + shape1 * (x/scale)^shape2) >= 0) & !is.na(x),
                 log(shape2) - shape2*log(scale) + (shape2-1)*log(x) + log1p(shape1 * (x / scale)^shape2),
                 NA))
  if (!log) {
    haz <- exp(haz)
  }
  return(haz)
}

#' Hazard function of the Weibull distribution
#'
#' @param x vector of quantiles
#' @param shape shape parameter
#' @param scale scale parameter, default to 1
#' @param log logical; if \code{TRUE}, returns the log hazard
#' @keywords internal
#' @return a vector of (log)-hazard
#' @export
hweibull <- function(x, shape, scale = 1, log = FALSE){
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "weibull")
  haz <- rep(NA, length(x))
  index <- which(x >=0 & !is.na(x))
  haz[index] <- -shape * log(scale) + log(shape) + (shape - 1) * log(x[index])
  if (!log) {
    haz <- exp(haz)
  }
  return(haz)
}

#' Quantile function of the extended Weibull distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pextweibull
#' @return vector of quantiles
#' @export
#' @keywords internal
#' @return a vector of quantiles
#' @rdname extweibull
qextweibull <- function(p,
                        scale = 1,
                        shape1 = 0,
                        shape2 = 1,
                        lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("\"p' must contain probabilities in (0,1)")
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extweibull")
  if (lower.tail) {
    p <- 1 - p
  }
  if (isTRUE(all.equal(c(shape1, shape2), c(0, 1), check.attributes = FALSE))) {
    return(-scale * log(p))
  } else if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))) {
    return(qweibull(
      p = p,
      shape = shape2,
      scale = scale,
      lower.tail = FALSE
    ))
  } else if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    return(scale * (p ^ (-shape1) - 1) / shape1)
  } else{
    return(scale * ((p ^ (-shape1) - 1) / shape1) ^ (1 / shape2))
  }
}



#' Distribution function of the Gompertz distribution
#'
#' @inheritParams pgpd
#' @export
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @keywords internal
#' @rdname gomp
pgomp <- function(q,
                  scale = 1,
                  shape = 0,
                  lower.tail = TRUE,
                  log.p = FALSE) {
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    # Exponential data
    return(pexp(
      q = q,
      rate = 1 / scale,
      lower.tail = lower.tail,
      log.p = log.p
    ))
  } else {
    p <-
      pmax(0, 1 - exp(-expm1(exp(
        log(shape) + log(pmax(0,q)) - log(scale)
      )) / shape))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Quantile function of the Gompertz distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @export
#' @return a vector of quantiles
#' @keywords internal
#' @rdname gomp
qgomp <- function(p,
                  scale = 1,
                  shape = 0,
                  lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("\"p' must contain probabilities in (0,1)")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    # Exponential data
    return(qexp(
      p = p,
      rate = 1 / scale,
      lower.tail = lower.tail
    ))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  return(scale / shape * log(1 - shape * log(1 - p)))
}

#' Distribution function of the Gompertz-Makeham distribution
#'
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @keywords internal
#' @rdname gompmake
pgompmake <- function(q,
                      scale = 1,
                      shape = 0,
                      lambda = 0,
                      lower.tail = TRUE,
                      log.p = FALSE) {
  check_elife_dist(scale = scale,
                   rate = lambda,
                   shape = shape,
                   family = "gompmake")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    # Exponential data
    return(pexp(
      q = q,
      rate = lambda + 1 / scale,
      lower.tail = lower.tail,
      log.p = log.p
    ))
  } else {
    p <- pmax(0, 1 - exp(-lambda * q - (exp(shape * q / scale) - 1) / shape))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Quantile function of the Gompertz-Makeham distribution
#'
#' @note The quantile function is defined in terms of Lambert's W function. Particular parameter combinations (small values of \code{lambda} lead to numerical overflow; the function throws a warning when this happens.
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return a vector of quantiles
#' @keywords internal
#' @rdname gompmake
qgompmake <- function(p,
                      scale = 1,
                      shape = 0,
                      lambda = 0,
                      lower.tail = TRUE) {
  check_elife_dist(scale = scale,
                   rate = lambda,
                   shape = shape,
                   family = "gompmake")
  # stopifnot("Install package \"gsl\" to use \"qgompmake\" with the Gompertz model.\n Try \"install.packages(\"gsl\")\"" = requireNamespace("gsl", quietly = TRUE))
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("\"p' must contain probabilities in (0,1)")
  }
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    # Exponential data - but parameter not identifiable
    return(qexp(
      p = p,
      rate = lambda + 1 / scale,
      lower.tail = lower.tail
    ))
  }
  if (isTRUE(all.equal(lambda, 0, check.attributes = FALSE))) {
    # Exponential data - but parameter not identifiable
    return(qgomp(
      p = p,
      scale = scale,
      shape = shape,
      lower.tail = lower.tail
    ))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  q <- rep(NA, length(p))
  y <-
    exp(1 / (scale * lambda)) * (1 - p) ^ (-shape / (scale * lambda)) / (scale *
                                                                           lambda)
  index <- is.finite(y)
  if (sum(index) != sum(is.finite(p) & p < 1)) {
    warning("Numerical overflow: some quantiles map to infinity.")
  }
  q[index] <-
    1 / (shape * lambda) - log(1 - p[index]) / lambda - scale / shape * .LambertW0(y[index])
  q[is.infinite(y)] <- Inf
  return(q)
}

#' Density function of the Gompertz-Makeham distribution
#'
#' @param x vector of quantiles.
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return a vector of density
#' @keywords internal
#' @rdname gompmake
dgompmake <- function(x,
                      scale = 1,
                      shape = 0,
                      lambda = 0,
                      log = FALSE) {
  check_elife_dist(scale = scale,
                   rate = lambda,
                   shape = shape,
                   family = "gompmake")
  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    # Exponential data - but parameter not identifiable
    return(dexp(
      x = x,
      rate = lambda + 1 / scale,
      log = log
    ))
  }
  ld1 <- log(lambda + exp(shape * x / scale) / scale)
  ld1 <- ifelse(is.finite(ld1), ld1, 0)
  ldens <- ld1 - lambda * x - (exp(shape * x / scale) - 1) / shape
  ldens <- ifelse(x < 0 | is.infinite(x), -Inf, ldens)
  if (log) {
    return(ldens)
  } else{
    return(exp(ldens))
  }
}

#' Distribution function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#'
pextgp <- function(q,
                   scale = 1,
                   shape1 = 0,
                   shape2 = 0,
                   lower.tail = TRUE,
                   log.p = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
  if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))) {
    # Exponential data
    return(pgpd(
      q = q,
      scale = scale,
      shape = shape2,
      lower.tail = lower.tail,
      log.p = log.p
    ))
  } else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    return(pgomp(
      q = q,
      scale = scale,
      shape = shape1,
      lower.tail = lower.tail,
      log.p = log.p
    ))
  }
  p <-
    1 - pmax(0, (1 + shape2 * (exp(
      shape1 * pmax(0, q) / scale
    ) - 1) / shape1)) ^ (-1 / shape2)
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}
#' Density function of the Gompertz distribution
#'
#' @param x vector of quantiles
#' @param scale positive scale parameter
#' @param shape non-negative shape parameter
#' @param log logical; if \code{TRUE}, return the log density
#' @export
#' @keywords internal
#' @return a vector of (log)-density.
#' @rdname gomp
dgomp <- function(x,
                  scale = 1,
                  shape = 0,
                  log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if (shape < 1e-8) {
    return(dexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  } else{
    ldens <-
      ifelse(x < 0 |
               is.infinite(x),
             -Inf,
             -log(scale) + (shape * x / scale - exp(shape * x / scale) / shape + 1 /
                              shape))
  }
  if (log) {
    return(ldens)
  } else{
    exp(ldens)
  }
}

#' Hazard function of the Gompertz distribution
#'
#' @inheritParams pgomp
#' @param log logical; if \code{TRUE}, return the log hazard
#' @export
#' @keywords internal
#' @rdname gomp
#' @return a vector of (log)-hazard.
hgomp <- function(x,
                  scale = 1,
                  shape = 0,
                  log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if (shape < 1e-8) {
    return(hexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  }
  haz <- ifelse(x < 0 | is.na(x),
             NA,
             -log(scale) + (shape * x / scale))
  if (!log) {
    haz <- exp(haz)
  }
  return(haz)
}

#' Hazard function of the Gompertz-Makeham distribution
#'
#' @inheritParams pgompmake
#' @param log logical; if \code{TRUE}, return the log hazard
#' @export
#' @keywords internal
#' @return a vector of (log)-hazard.
#' @rdname gompmake
hgompmake <- function(x,
                  scale = 1,
                  shape = 0,
                  lambda = 0,
                  log = FALSE) {
  check_elife_dist(rate = lambda,
                   scale = scale,
                   shape = shape,
                   family = "gompmake")
  if(lambda == 0){
    return(hgomp(x = x, scale = scale, shape = shape, log = log))
  }
  if (shape < 1e-8) {
    return(hexp(
      x = x,
      rate = lambda + 1 / scale,
      log = log
    ))
  }
  haz <- ifelse(x < 0 | is.na(x),
                NA,
                lambda + exp(shape * x / scale) / scale)
  if (log) {
    haz <- log(haz)
  }
  return(haz)
}

#' Density function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @param log logical; if \code{TRUE}, return the log-density
#' @return a vector of (log)-density of the same length as \code{x}
#' @export
#' @keywords internal
#' @rdname extgp
dextgp <- function(x,
                   scale = 1,
                   shape1 = 0,
                   shape2 = 0,
                   log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
  if (abs(shape2) < 1e-8 && abs(shape1) < 1e-8) {
    return(dexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  } else if (abs(shape2) < 1e-8 && abs(shape1) > 1e-8) {
    #Gompertz
    ldens <-
      -log(scale) + (shape1 * x / scale - exp(shape1 * x / scale) / shape1 + 1 /
                       shape1)
    ldens <- ifelse(x < 0, -Inf, ldens)
  } else if (abs(shape2) >= 1e-8 &&
             abs(shape1) < 1e-8) {
    #generalized Pareto
    return(dgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape2,
      log = log
    ))
  } else {
    #extended
    ldens <-
      suppressWarnings(-log(scale) + (-1 / shape2 - 1) * log(shape2 * (exp(shape1 *
                                                                             x / scale) - 1) / shape1 + 1) + shape1 * x / scale)

    ldens <- ifelse(x < 0 | is.infinite(x), -Inf, ldens)
    if (shape2 < 0) {
      ldens <- ifelse(x > scale * log(1 - shape1 / shape2) / shape1,
                      -Inf,
                      ldens)
    }
  }
  if (log) {
    return(ldens)
  } else{
    exp(ldens)
  }
}

#' Hazard function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @param log logical; if \code{TRUE}, return the log hazard
#' @return a vector of (log)-hazard of the same length as \code{x}
#' @export
#' @keywords internal
#' @rdname extgp
hextgp <- function(x,
                   scale = 1,
                   shape1 = 0,
                   shape2 = 0,
                   log = FALSE) {
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
  if (abs(shape2) < 1e-8 && abs(shape1) < 1e-8) {
    return(hexp(
      x = x,
      rate = 1 / scale,
      log = log
    ))
  } else if (abs(shape2) < 1e-8 && abs(shape1) > 1e-8) {
    #Gompertz
    return(hgomp(
      x = x,
      scale = scale,
      shape = shape1,
      log = log
    ))
  } else if (abs(shape2) >= 1e-8 &&
             abs(shape1) < 1e-8) {
    #generalized Pareto
    return(hgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape2,
      log = log
    ))
  } else{
    # extended generalized Pareto
    haz <- suppressWarnings(
      ifelse(x < 0 | is.na(x) & (shape1 + shape2 * (exp(shape1 * x / scale) - 1)) < 0,
                  NA,
                  log(shape1) - log(scale) + shape1 * x / scale - log(shape1 + shape2 * (exp(shape1 * x / scale) - 1)))
                  )
  if (!log) {
    haz <- exp(haz)
  }
  return(haz)
  }
}



#' Quantile function of the extended generalized Pareto distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @inheritParams pextgp
#' @return a vector of quantiles
#' @export
#' @keywords internal
#' @rdname extgp
qextgp <- function(p,
                   scale = 1,
                   shape1 = 0,
                   shape2 = 0,
                   lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("\"p\" must contain probabilities in [0,1]")
  }
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
  if (isTRUE(all.equal(shape1, 0, check.attributes = FALSE))) {
    # Exponential data
    return(qgpd(
      p = p,
      scale = scale,
      shape = shape2,
      lower.tail = lower.tail
    ))
  } else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    return(qgomp(
      p = p,
      scale = scale,
      shape = shape1,
      lower.tail = lower.tail
    ))
  }
  if (lower.tail) {
    p <- 1 - p
  }
  return(scale / shape1 * log(shape1 / shape2 * (p ^ (-shape2) - 1) + 1))
}

#' Distribution function of the Perks-Makeham distribution
#'
#' @param q vector of quantiles.
#' @param rate rate parameter (\eqn{\nu})
#' @param shape shape parameter (\eqn{\alpha})
#' @param lambda exponential rate of the Makeham component
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname perksmake
pperksmake <- function(q,
                       rate = 1,
                       shape = 1,
                       lambda = 0,
                       lower.tail = TRUE,
                       log.p = FALSE) {
  q <- pmax(q, 0)
  check_elife_dist(rate = c(rate, lambda),
                   shape = shape,
                   family = "perksmake")
  if (rate == 0) {
    p <- 1 - exp(-(lambda + shape / (shape + 1)) * q)
  } else {
    p <-
      1 - exp(-lambda * q + (log1p(shape) - log1p(shape * exp(rate * q))) /
                rate)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Density function of the Perks-Makeham distribution
#'
#' @inheritParams pperksmake
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname perksmake
dperksmake <- function (x,
                        rate = 1,
                        shape = 1,
                        lambda = 0,
                        log = FALSE) {
  check_elife_dist(rate = c(rate, lambda),
                   shape = shape,
                   family = "perksmake")
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    return(dexp(
      x = x,
      rate = lambda + shape / (shape + 1),
      log = log
    ))
  }
  d <-  log1p(shape) / rate +
    log(lambda + shape * exp(rate * x) /
              (1 + shape *exp(rate * x))) - lambda * x - log1p(shape * exp(rate * x)) / rate
  d[!is.finite(x) | x < 0] <- -Inf
  d[is.na(x)] <- NA
  if (!log) {
    d <- exp(d)
  }
  return(d)
}


#' Hazard function of the Perks-Makeham distribution
#'
#' @inheritParams pperksmake
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname perksmake
hperksmake <- function (x,
                        rate = 1,
                        shape = 1,
                        lambda = 0,
                        log = FALSE) {
  check_elife_dist(rate = c(rate, lambda),
                   shape = shape,
                   family = "perksmake")
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    return(hexp(
      x = x,
      rate = lambda + shape / (shape + 1),
      log = log
    ))
  }
  haz <- ifelse(x < 0 | is.na(x),
                NA,
                lambda + shape * exp(rate*x) / (1 + shape * exp(rate*x)))
  if (log) {
    haz <- log(haz)
  }
  return(haz)
}


#' Hazard function of the Perks distribution
#'
#' @inheritParams pperks
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname perks
hperks <- function (x,
                    rate = 1,
                    shape = 1,
                    log = FALSE) {
  check_elife_dist(rate = c(rate, 0),
                   shape = shape,
                   family = "perksmake")
  return(hperksmake(x = x, rate = rate, shape = shape, log = log, lambda = 0))
}

#' Quantile function of the Perks-Makeham distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pperksmake
#' @return vector of quantiles
#' @export
#' @importFrom stats uniroot
#' @return a vector of quantiles
#' @keywords internal
#' @rdname perksmake
qperksmake <- function(p,
                       rate = 1,
                       shape = 1,
                       lambda = 0,
                       lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("\"p' must contain probabilities in (0,1)")
  check_elife_dist(rate = c(rate, lambda),
                   shape = shape,
                   family = "perksmake")
  if (lower.tail) {
    p <- 1 - p
  }
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    rate_e <- lambda + shape / (1 + shape)
    return(qexp(p, rate = rate_e, lower.tail = FALSE))
  } else if (isTRUE(all.equal(lambda, 0, check.attributes = FALSE))) {
    return(log((p ^ (-rate) * (1 + shape) - 1) / shape) / rate)
  } else{
    # No closed-form expression, resort to numerical root finding
    return(sapply(
      log(p) - log1p(shape) / rate,
      FUN = function(q) {
        suppressWarnings(uniroot(
          f = function(x, q) {
            -lambda * x - log1p(shape * exp(rate * x)) / rate - q
          },
          q = q,
          interval = c(0, 1e20)
        )$root)
      }
    ))
  }
}

#' Distribution function of the Perks distribution
#'
#' @param q vector of quantiles.
#' @param rate rate parameter (\eqn{\nu})
#' @param shape shape parameter (\eqn{\alpha})
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname perks
pperks <- function(q,
                   rate = 1,
                   shape = 1,
                   lower.tail = TRUE,
                   log.p = FALSE) {
  pperksmake(
    q = q,
    rate = rate,
    shape = shape,
    lambda = 0,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' Density function of the Perks distribution
#'
#' @inheritParams pperks
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname perks
dperks <- function (x,
                    rate = 1,
                    shape = 1,
                    log = FALSE) {
  dperksmake(
    x = x,
    rate = rate,
    shape = shape,
    log = log
  )
}

#' Quantile function of the Perks distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pperks
#' @return a vector of quantiles
#' @export
#' @keywords internal
#' @rdname perks
qperks <- function(p,
                   rate = 1,
                   shape = 1,
                   lower.tail = TRUE) {
  qperksmake(
    p = p,
    rate = rate,
    shape = shape,
    lower.tail = lower.tail
  )
}

#' Beard distribution
#'
#' Distribution function, density, hazard, quantile function and random number generation for the Beard distribution
#'
#' @param n sample size
#' @param rate rate parameter (\eqn{\nu})
#' @param shape1 shape parameter (\eqn{\alpha})
#' @param shape2 shape parameter (\eqn{\beta})
#' @return a vector of random variates
#' @export
#' @keywords internal
#' @rdname beard
rbeard <- function(n, rate, shape1, shape2){
  n <- as.integer(n)
  stopifnot(n >= 1)
  qbeard(runif(n), rate = rate, shape1 = shape1, shape2 = shape2)
}


#' @param q vector of quantiles.
#' @param rate rate parameter (\eqn{\nu})
#' @param shape1 shape parameter (\eqn{\alpha})
#' @param shape2 shape parameter (\eqn{\beta})
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname beard
pbeard <- function(q,
                   rate = 1,
                   shape1 = 1,
                   shape2 = 1,
                   lower.tail = TRUE,
                   log.p = FALSE) {
  pbeardmake(
    q = q,
    rate = rate,
    shape1 = shape1,
    shape2 = shape2,
    lambda = 0,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @inheritParams pbeard
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname beard
dbeard <- function (x,
                    rate = 1,
                    shape1 = 1,
                    shape2 = 1,
                    log = FALSE) {
  dbeardmake(
    x = x,
    rate = rate,
    shape1 = shape1,
    shape2 = shape2,
    lambda = 0,
    log = log
  )
}

#' @inheritParams pbeard
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname beard
hbeard <- function (x,
                    rate = 1,
                    shape1 = 1,
                    shape2 = 1,
                    log = FALSE) {
  hbeardmake(
    x = x,
    rate = rate,
    shape1 = shape1,
    shape2 = shape2,
    lambda = 0,
    log = log
  )
}

#' @param p vector of probabilities.
#' @inheritParams pbeard
#' @return a vector of quantiles
#' @export
#' @keywords internal
#' @rdname beard
qbeard <- function(p,
                   rate = 1,
                   shape1 = 1,
                   shape2 = 1,
                   lower.tail = TRUE) {
  qbeardmake(
    p = p,
    rate = rate,
    shape1 = shape1,
    shape2 = shape2,
    lambda = 0,
    lower.tail = lower.tail
  )
}

#' Distribution function of the Beard-Makeham distribution
#'
#' @param q vector of quantiles.
#' @param rate shape parameter (\eqn{\nu})
#' @param shape1 shape parameter (\eqn{\alpha})
#' @param shape2 shape parameter (\eqn{\beta})
#' @param lambda exponential rate of the Makeham component
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
#' @rdname beardmake
pbeardmake <- function(q,
                       rate = 1,
                       shape1 = 1,
                       shape2 = 1,
                       lambda = 0,
                       lower.tail = TRUE,
                       log.p = FALSE) {
  q <- pmax(q,0)
  if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(
      pperksmake(
        q = q,
        rate = rate,
        shape = shape1,
        lambda = lambda,
        lower.tail = lower.tail,
        log.p = log.p
      )
    )
  } else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(
      pgompmake(
        q = q,
        scale = 1 / shape1,
        shape = rate / shape1,
        lambda = lambda,
        lower.tail = lower.tail,
        log.p = log.p
      )
    )
  }
  check_elife_dist(
    rate = c(rate, lambda),
    shape = c(shape1, shape2),
    family = "beardmake"
  )
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    p <- 1 - exp(-(lambda + shape1 / (shape1 * shape2 + 1)) * q)
  } else {
    p <-
      1 - exp((log1p(shape1 * shape2) - log1p(shape1 * shape2 * exp(rate *q))) / (shape2 * rate) - lambda * q)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' Density function of the Beard-Makeham distribution
#'
#' @inheritParams pbeardmake
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
#' @rdname beardmake
dbeardmake <- function (x,
                        rate = rate,
                        shape1 = 1,
                        shape2 = 1,
                        lambda = 0,
                        log = FALSE) {
  if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(dperksmake(
      x = x,
      rate = rate,
      shape = shape1,
      lambda = lambda,
      log = log
    ))
  } else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(dgompmake(
      x = x,
        scale = 1 / shape1,
        shape = rate / shape1,
        lambda = lambda,
      log = log
    ))
  }
  check_elife_dist(
    rate = c(rate, lambda),
    shape = c(shape1, shape2),
    family = "beardmake"
  )
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    return(dexp(
      x = x,
      rate = lambda + shape1 / (1 + shape1 * shape2),
      log = log
    ))
  }
  d <-
    (log1p(shape1 * shape2) - log1p(shape1 * shape2 * exp(rate * x))) /
    (shape2 * rate) - lambda * x + log(lambda + shape1 * exp(rate * x) / (1 + shape1 *
                                                                                                                                               shape2 * exp(rate * x)))
  d[!is.finite(x) | x < 0] <- -Inf
  d[is.na(x)] <- NA
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

#' Hazard function of the Beard-Makeham distribution
#'
#' @inheritParams pbeardmake
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the hazard
#' @return a vector of (log)-hazard.
#' @export
#' @keywords internal
#' @rdname beardmake
hbeardmake <- function (x,
                        rate = rate,
                        shape1 = 1,
                        shape2 = 1,
                        lambda = 0,
                        log = FALSE) {
  if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(hperksmake(
      x = x,
      rate = rate,
      shape = shape1,
      lambda = lambda,
      log = log
    ))
  } else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(hgompmake(
      x = x,
      scale = 1 / shape1,
      shape = rate / shape1,
      lambda = lambda,
      log = log
    ))
  }
  check_elife_dist(
    rate = c(rate, lambda),
    shape = c(shape1, shape2),
    family = "beardmake"
  )
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    return(hexp(
      x = x,
      rate = lambda + shape1 / (1 + shape1 * shape2),
      log = log
    ))
  }
  haz <- ifelse(x < 0 | is.na(x),
                NA,
                lambda + shape1 * exp(rate*x) / (1 + shape1 * shape2 * exp(rate*x)))
  if (log) {
    haz <- log(haz)
  }
  return(haz)
}

#' Quantile function of the Beard-Makeham distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pbeardmake
#' @return a vector of quantiles
#' @export
#' @keywords internal
#' @rdname beardmake
qbeardmake <- function(p,
                       rate = 1,
                       shape1 = 1,
                       shape2 = 1,
                       lambda = 0,
                       lower.tail = TRUE) {
  if (isTRUE(all.equal(shape2, 1, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(
      qperksmake(
        p = p,
        rate = rate,
        shape = shape1,
        lambda = lambda,
        lower.tail = lower.tail
      )
    )
  }  else if (isTRUE(all.equal(shape2, 0, check.attributes = FALSE))) {
    # If subcase, return this instead
    return(
      qgompmake(
        p = p,
        scale = 1 / shape1,
        shape = rate / shape1,
        lambda = lambda,
        lower.tail = lower.tail
      )
    )
  }
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("\"p' must contain probabilities in (0,1)")
  check_elife_dist(
    rate = c(rate, lambda),
    shape = c(shape1, shape2),
    family = "beardmake"
  )
  if (lower.tail) {
    p <- 1 - p
  }
  if (isTRUE(all.equal(rate, 0, check.attributes = FALSE))) {
    rate_e <- lambda + shape1 / (1 + shape1 * shape2)
    return(qexp(p, rate = rate_e, lower.tail = FALSE))
  } else if (isTRUE(all.equal(lambda, 0, check.attributes = FALSE))) {
    return(log((
      p ^ (-rate * shape2) * (1 + shape1 * shape2) - 1
    ) / (shape1 * shape2)) / rate)
  } else{
    # No closed-form expression, resort to numerical root finding
    return(sapply(
      log(p) - log1p(shape2 * shape1) / (rate * shape2),
      FUN = function(q) {
        suppressWarnings(uniroot(
          f = function(x, q) {
            -lambda * x - log1p(shape2 * shape1 * exp(rate * x)) / (rate * shape2) - q
          },
          q = q,
          interval = c(0, 1e20)
        )$root)
      }
    ))
  }
}


#' Excess lifetime distributions
#'
#' Quantile, distribution, density and hazard functions of excess lifetime distribution
#' for threshold exceedances.
#' @param q vector of quantiles.
#' @param p vector of probabilities
#' @param n sample size
#' @param rate rate parameter(s); for models with Makeham component, the last entry should be part of the rate vector
#' @param scale scale parameter
#' @param shape vector of shape parameter(s).
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log,log.p logical; if \code{TRUE}, values are returned on the logarithmic scale (default to \code{FALSE}).
#' @param family string indicating the parametric model, one of \code{exp}, \code{gp}, \code{gomp}, \code{gompmake}, \code{weibull}, \code{extgp}, \code{extweibull}, \code{perks}, \code{perksmake}, \code{beard} and \code{beardmake}
#' @name elife
#' @returns depending on the function type, a vector of probabilities (\code{pelife}), survival probabilities (\code{selife}), quantiles (\code{qelife}), density (\code{delife}), or hazard (\code{helife}). The function \code{relife} returns a random sample of size \code{n} from the distribution.
NULL

#' @rdname elife
#' @export
#' @keywords internal
qelife <- function(p,
                   rate,
                   scale,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   ),
                   lower.tail = TRUE) {
  if(family == "gpd"){ family <- "gp"}
  family <- match.arg(family)
  if (missing(shape) & family != "exp") {
    stop("Missing \"shape\" parameter.")
  }
  if(family == "exp"){
    if(missing(rate) & missing(scale)){
      stop("No scale parameter is provided")
    } else if(missing(scale)){
      scale <- 1/rate
    }
  }
  # Technically, this is called again in the other functions...
  if(family %in%  c("exp","weibull")){
  check_elife_dist(rate = rate,
                   scale = scale,
                   shape = shape,
                   family = family)
  }
  switch(
    family,
    exp = qexp(p, rate = 1 / scale, lower.tail = lower.tail),
    gp = qgpd(
      p = p,
      loc = 0,
      scale = scale,
      shape = shape,
      lower.tail = lower.tail
    ),
    weibull = qweibull(
      p = p,
      shape = shape,
      scale = scale,
      lower.tail = lower.tail
    ),
    gomp = qgomp(
      p = p,
      scale = scale,
      shape = shape,
      lower.tail = lower.tail
    ),
    gompmake = qgompmake(
      p = p,
      scale = scale,
      lambda = rate,
      shape = shape,
      lower.tail = lower.tail
    ),
    extgp = qextgp(
      p = p,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      lower.tail = lower.tail
    ),
    extweibull = qextweibull(
      p = p,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      lower.tail = lower.tail
    ),
    perks = qperks(
      p = p,
      rate = rate,
      shape = shape,
      lower.tail = lower.tail
    ),
    perksmake = qperksmake(
      p = p,
      lambda = rate[2],
      shape = shape,
      rate = rate[1],
      lower.tail = lower.tail
    ),
    beard = qbeard(
      p = p,
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate,
      lower.tail = lower.tail
    ),
    beardmake = qbeardmake(
      p = p,
      lambda = rate[2],
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate[1],
      lower.tail = lower.tail
    )
  )
}

#' @rdname elife
#' @export
pelife <- function(q,
                   rate,
                   scale,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   ),
                   lower.tail = TRUE,
                   log.p = FALSE) {
  if(family == "gpd"){ family <- "gp"}
  family <- match.arg(family)
  if (missing(shape) & family != "exp") {
    stop("Missing \"shape\" parameter.")
  }
  if(family == "exp"){
    if(missing(rate) & missing(scale)){
      stop("No scale parameter is provided")
    } else if(missing(scale)){
      scale <- 1/rate
    }
  }
  if(family %in%  c("exp","weibull")){
    check_elife_dist(rate = rate,
                     scale = scale,
                     shape = shape,
                     family = family)
  }
  switch(
    family,
    exp = pexp(
      q = q,
      rate = 1 / scale,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    gp = pgpd(
      q = q,
      loc = 0,
      scale = scale,
      shape = shape,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    weibull = pweibull(
      q = q,
      shape = shape,
      scale = scale,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    gomp = pgomp(
      q = q,
      scale = scale,
      shape = shape,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    gompmake = pgompmake(
      q = q,
      scale = scale,
      shape = shape,
      lambda = rate,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    extgp = pextgp(
      q = q,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      lower.tail = lower.tail,
      log.p = log.p
    ),
    extweibull = pextweibull(
      q = q,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      lower.tail = lower.tail,
      log.p = log.p
    ),
    perks = pperks(
      q = q,
      rate = rate,
      shape = shape,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    perksmake = pperksmake(
      q = q,
      rate = rate[1],
      lambda = rate[2],
      shape = shape,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    beard = pbeard(
      q = q,
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate,
      lower.tail = lower.tail,
      log.p = log.p
    ),
    beardmake = pbeardmake(
      q = q,
      rate = rate[1],
      lambda = rate[2],
      shape1 = shape[1],
      shape2 = shape[2],
      lower.tail = lower.tail,
      log.p = log.p
    )
  )
}

#' @rdname elife
#' @export
selife <- function(q,
                   rate,
                   scale,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   ),
                   log.p = FALSE) {
  if(family == "gpd"){ family <- "gp"}
  pelife(q = q,
         rate = rate,
         scale = scale,
         shape = shape,
         family = family,
         lower.tail = FALSE,
         log.p = log.p)
}

#' @rdname elife
#' @export
relife <- function(n,
                   scale = 1,
                   rate,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   )) {
  if(family == "gpd"){ family <- "gp"}
  family <- match.arg(family)
  if (missing(shape) & family != "exp") {
    stop("Missing \"shape\" parameter.")
  }
  if(family == "exp"){
    if(missing(rate) & missing(scale)){
      stop("No scale parameter is provided")
    } else if(!missing(rate) & missing(scale)){
      scale <- 1/rate
    }
  }
  if(family %in%  c("exp","weibull")){
    check_elife_dist(rate = rate,
                     scale = scale,
                     shape = shape,
                     family = family)
  }
  switch(
    family,
    exp = rexp(n = n, rate = 1 / scale),
    gp = qgpd(
      p = runif(n),
      loc = 0,
      scale = scale,
      shape = shape
    ),
    weibull = rweibull(n = n, shape = shape, scale = scale),
    gomp = qgomp(
      p = runif(n),
      scale = scale,
      shape = shape
    ),
    gompmake = qgompmake(
      p = runif(n),
      scale = scale,
      lambda = rate,
      shape = shape
    ),
    extgp = qextgp(
      p = runif(n),
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2]
    ),
    extweibull = qextweibull(
      p = runif(n),
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2]
    ),
    perks = qperks(
      p = runif(n),
      shape = shape,
      rate = rate
    ),
    perksmake = qperksmake(
      p = runif(n),
      lambda = rate[2],
      shape = shape,
      rate = rate[1]
    ),
    beard = qbeard(
      p = runif(n),
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate
    ),
    beardmake = qbeardmake(
      p = runif(n),
      lambda = rate[2],
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate[1]
    )
  )
}


#' @rdname elife
#' @export
delife <- function(x,
                   scale = 1,
                   rate,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   ),
                   log = FALSE) {
  if(family == "gpd"){ family <- "gp"}
  family <- match.arg(family)
  if (missing(shape) & family != "exp") {
    stop("Missing \"shape\" parameter.")
  }
  if(family == "exp"){
    if(missing(rate) & missing(scale)){
      stop("No scale parameter is provided")
    } else if(!missing(rate) & missing(scale)){
      scale <- 1/rate
    }
  }
  if(family %in%  c("exp","weibull")){
    check_elife_dist(rate = rate,
                     scale = scale,
                     shape = shape,
                     family = family)
  }
  switch(
    family,
    exp = dexp(
      x = x,
      rate = 1 / scale,
      log = log
    ),
    gp = dgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape,
      log = log
    ),
    weibull = dweibull(
      x = x,
      shape = shape,
      scale = scale,
      log = log
    ),
    gomp = dgomp(
      x = x,
      scale = scale,
      shape = shape,
      log = log
    ),
    gompmake = dgompmake(
      x = x,
      scale = scale,
      lambda = rate,
      shape = shape,
      log = log
    ),
    extgp = dextgp(
      x = x,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      log = log
    ),
    extweibull = dextweibull(
      x = x,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      log = log
    ),
    perks = dperks(
      x = x,
      shape = shape,
      rate = rate,
      log = log
    ),
    perksmake = dperksmake(
      x = x,
      lambda = rate[2],
      shape = shape,
      rate = rate[1],
      log = log
    ),
    beard = dbeard(
      x = x,
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate,
      log = log
    ),
    beardmake =  dbeardmake(
      x = x,
      lambda = rate[2],
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate[1],
      log = log
    ),
  )
}


#' @rdname elife
#' @export
helife <- function(x,
                   scale = 1,
                   rate,
                   shape,
                   family = c(
                     "exp",
                     "gp",
                     "weibull",
                     "gomp",
                     "gompmake",
                     "extgp",
                     "extweibull",
                     "perks",
                     "perksmake",
                     "beard",
                     "beardmake"
                   ),
                   log = FALSE) {
  if(family == "gpd"){ family <- "gp"}
  family <- match.arg(family)
  if (missing(shape) & family != "exp") {
    stop("Missing \"shape\" parameter.")
  }
  if(family == "exp"){
    if(missing(rate) & missing(scale)){
      stop("No scale parameter is provided")
    } else if(!missing(rate) & missing(scale)){
      scale <- 1/rate
    }
  }
  if(family %in%  c("exp","weibull")){
    check_elife_dist(rate = rate,
                     scale = scale,
                     shape = shape,
                     family = family)
  }
  switch(
    family,
    exp = hexp(
      x = x,
      rate = 1 / scale,
      log = log
    ),
    gp = hgpd(
      x = x,
      loc = 0,
      scale = scale,
      shape = shape,
      log = log
    ),
    weibull = hweibull(
      x = x,
      shape = shape,
      scale = scale,
      log = log
    ),
    gomp = hgomp(
      x = x,
      scale = scale,
      shape = shape,
      log = log
    ),
    gompmake = hgompmake(
      x = x,
      scale = scale,
      lambda = rate,
      shape = shape,
      log = log
    ),
    extgp = hextgp(
      x = x,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      log = log
    ),
    extweibull = hextweibull(
      x = x,
      scale = scale,
      shape1 = shape[1],
      shape2 = shape[2],
      log = log
    ),
    perks = hperks(
      x = x,
      shape = shape,
      rate = rate,
      log = log
    ),
    perksmake = hperksmake(
      x = x,
      lambda = rate[2],
      shape = shape,
      rate = rate[1],
      log = log
    ),
    beard = hbeard(
      x = x,
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate,
      log = log
    ),
    beardmake =  hbeardmake(
      x = x,
      lambda = rate[2],
      shape1 = shape[1],
      shape2 = shape[2],
      rate = rate[1],
      log = log
    ),
  )
}


#' Check parameters of extended lifetime distributions
#' @inheritParams pelife
#' @export
#' @return The function has no return value and is only used to throw error for invalid arguments.
#' @keywords internal
check_elife_dist <- function(rate,
                             scale,
                             shape,
                             family = c(
                               "exp",
                               "gp",
                               "weibull",
                               "gomp",
                               "gompmake",
                               "extgp",
                               "extweibull",
                               "perks",
                               "beard",
                               "perksmake",
                               "beardmake"
                             )) {
  family <- match.arg(family)
  if (family == "gp") {
    stopifnot("Invalid scale parameter: must be a positive scalar." =
                isTRUE(all(
                  length(scale) == 1L,
                  is.finite(scale),
                  scale > 0
                )))
    stopifnot(
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be finite." = isTRUE(is.finite(shape))
    )
  } else if (family == "weibull") {
    stopifnot("Invalid scale parameter: must be a positive scalar." =
                isTRUE(all(
                  length(scale) == 1L,
                  is.finite(scale),
                  scale > 0
                )))
    stopifnot(
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be positive." = isTRUE(is.finite(shape) &
                                               shape > 0)
    )
  } else if (family == "gomp") {
    stopifnot("Invalid scale parameter: must be a positive scalar." =
                isTRUE(all(
                  length(scale) == 1L,
                  is.finite(scale),
                  scale > 0
                )))
    stopifnot(
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be non-negative." = isTRUE(is.finite(shape) &
                                                   shape >= 0)
    )
  } else if (family == "gompmake") {
      stopifnot(
      "Invalid scale parameter." =
        isTRUE(all(
          length(scale) == 1L,
          is.finite(scale),
          scale > 0
        )),
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be non-negative." = isTRUE(all(is.finite(shape), shape >= 0)),
      "\"rate\" must be non-negative." = isTRUE(all(is.finite(rate), rate >= 0))
    )
  } else if (family == "extgp") {
    stopifnot("Invalid scale parameter: must be a positive scalar." =
                isTRUE(all(
                  length(scale) == 1L,
                  is.finite(scale),
                  scale > 0
                )))
    stopifnot(
      "\"shape\" should be a vector of length 2." = length(shape) == 2L,
      "\"shape1\" must be non-negative." = isTRUE(is.finite(shape[1]) &
                                                    shape[1] >= 0),
      "\"shape2\" must be finite." = isTRUE(is.finite(shape[2]))
    )
  } else if (family == "extweibull") {
    stopifnot("Invalid scale parameter: must be a positive scalar." =
                isTRUE(all(
                  length(scale) == 1L,
                  is.finite(scale),
                  scale > 0
                )))
    stopifnot(
      "\"shape\" should be a vector of length 2." = length(shape) == 2L,
      "\"shape\" must be finite." = isTRUE(all(is.finite(shape))),
      "\"shape2\" must be positive." = shape[2] > 0
    )
  }  else if (family == "perks") {
    stopifnot(
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be finite." = isTRUE(is.finite(shape)),
      "\"shape\" must be positive." = isTRUE(shape > 0),
      "\"rate\" should be a vector of length 1." = length(rate) == 1L,
      "\"rate\" must be finite." = isTRUE(is.finite(rate)),
      "\"rate\" must be non-negative." = isTRUE(rate > 0)
    )
  }  else if (family == "perksmake") {
    stopifnot(
      "\"shape\" should be a vector of length 1." = length(shape) == 1L,
      "\"shape\" must be finite." = isTRUE(all(is.finite(shape))),
      "\"shape\" must be positive." = isTRUE(shape > 0)
    )
    stopifnot("Invalid rate parameter: must be a vector of length 2 of non-negative entries." =
                isTRUE(all(
                  length(rate) == 2L,
                  is.finite(rate),
                  isTRUE(all(rate >= 0))
                )))
  }   else if (family == "beard") {
    stopifnot(
      "\"shape\" should be a vector of length 2." = length(shape) == 2L,
      "\"shape\" must be finite." = isTRUE(all(is.finite(shape))),
      "First argument of \"shape\" must be positive." = isTRUE(shape[1] > 0),
      "Second argument of \"shape\" must be non-negative." = isTRUE(shape[2] > 0),
      "\"rate\" should be a vector of length 1." = length(rate) == 1L,
      "\"rate\" must be finite." = isTRUE(is.finite(rate)),
      "\"rate\" must be non-negative." = isTRUE(rate > 0)
    )
  }  else if (family == "beardmake") {
    stopifnot(
      "\"shape\" should be a vector of length 2." = length(shape) == 2L,
      "\"shape\" must be finite." = isTRUE(all(is.finite(shape))),
      "First argument of \"shape\" must be positive." = isTRUE(shape[1] > 0),
      "Second argument of \"shape\" must be non-negative." = isTRUE(shape[2] > 0)
    )
    stopifnot("Invalid rate parameter: must be a vector of length 2 of non-negative entries." =
                isTRUE(all(
                  length(rate) == 2L,
                  isTRUE(all(is.finite(rate))),
                  isTRUE(all(rate >= 0))
                )))
  }
}


# The Lambert W function
#
# This is extracted from VGAM
# @author Thomas W. Yee
# @keywords internal
.LambertW0 <- function (x,
                        tolerance = 1e-10,
                        maxit = 50)
{
  ans <- x
  ans[!is.na(x) & x < -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >=
                                              -exp(-1)])
  ans[!is.na(x) & x >= 0] <- sqrt(x[!is.na(x) & x >= 0]) / 2
  cutpt <- 3
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])
    L2 <- log(L1)
    wzinit <- L1 - L2 + (L2 + (L2 * (-2 + L2) / (2) + (
      L2 *
        (6 + L2 * (-9 + L2 * 2)) /
        (6) + L2 * (-12 + L2 * (36 +
                                  L2 * (-22 + L2 * 3))) /
        (12 * L1)
    ) / L1) / L1) / L1
    ans[myTF] <- wzinit
  }
  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x) / (exp2 + exp1 - ((ans + 2) * (exp2 -
                                                         x) / (2 * (ans + 1))))
    ans <- ans - delta
    if (all(is.na(delta)) || max(abs(delta), na.rm = TRUE) <
        tolerance)
      break
    if (ii == maxit)
      warning("did not converge")
  }
  ans[is.infinite(NA)] <- Inf
  ans
}
