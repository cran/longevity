#' Likelihood for doubly interval truncated data
#'
#' Computes the log-likelihood for various parametric models suitable for threshold
#' exceedances. If threshold is non-zero, then only right-censored, observed event time and interval censored
#' data whose timing exceeds the thresholds are kept.
#' @keywords internal
#' @param par vector of parameters
#' @param thresh vector of thresholds
#' @inheritParams npsurv
#' @param ltrunc1 lower truncation limit, default to \code{NULL}
#' @param rtrunc1 upper truncation limit, default to \code{NULL}
#' @param ltrunc2 lower truncation limit, default to \code{NULL}
#' @param rtrunc2 upper truncation limit, default to \code{NULL}
#' @param weights weights for observations
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}), Gompertz-Makeham (\code{gompmake}) or extended generalized Pareto (\code{extgp}).
#' @param arguments a named list specifying default arguments of the function that are common to all \code{elife} calls
#' @param ... additional arguments for optimization, currently ignored.
#' @return log-likelihood value
#' @export
nll_ditrunc_elife <-
  function(par,
          time,
          ltrunc1 = NULL,
          rtrunc1 = NULL,
          ltrunc2 = NULL,
          rtrunc2 = NULL,
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
          thresh = 0,
          weights = rep(1, length(time)),
          arguments = NULL,
          ...){
    if(!is.null(arguments)){
      call <- match.call(expand.dots = FALSE)
      arguments <- check_arguments(func = nll_ditrunc_elife, call = call, arguments = arguments)
      return(do.call(nll_ditrunc_elife, args = arguments))
    }

  family <- match.arg(family)
  stopifnot("Threshold must be non-negative" = thresh >= 0,
            "Only a single threshold value is allowed" = family == "gppiece" | length(thresh) == 1L)
  # Format time entry and check for input
  par <- par[is.finite(par)]
  stopifnot("Incorrect `ltrunc1` argument." = is.null(ltrunc1) || (!is.null(ltrunc1) & length(ltrunc1) %in% c(1L, length(time))),
            "Incorrect `rtrunc1` argument." = is.null(rtrunc1) || (!is.null(rtrunc1) & length(rtrunc1) %in% c(1L, length(time))),
            "Incorrect `ltrunc2` argument." = is.null(ltrunc2) || (!is.null(ltrunc2) & length(ltrunc2) %in% c(1L, length(time))),
            "Incorrect `rtrunc2` argument." = is.null(rtrunc2) || (!is.null(rtrunc2) & length(rtrunc2) %in% c(1L, length(time))))
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- time > thresh[1]
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    if(!is.null(ltrunc1)){ #both ltrc and ltrt
      ltrunc1 <- pmax(0, ltrunc1[ind] - thresh[1])
    }
    if(!is.null(ltrunc2)){
      ltrunc2 <- pmax(0, ltrunc2[ind] - thresh[1])
    }
    if(!is.null(rtrunc1)){
      rtrunc1 <- pmax(0, rtrunc1[ind] - thresh[1])
    }
    if(!is.null(rtrunc2)){
      rtrunc2 <- pmax(0, rtrunc2[ind] - thresh[1])
    }
    # Sometimes, the thresholding leads to ltrunc1=rtrunc1=0
    bad <- which(rtrunc1 == 0)
    if(length(bad) > 0){
      ltrunc1[bad] <- ltrunc2[bad]
      rtrunc1[bad] <- rtrunc2[bad]
      ltrunc2[bad] <- NA
      rtrunc2[bad] <- NA
    }
  }
  # For gppiece model
  if(family == "gomp"){
    family <- "extgp"
    par <- c(par, 0)
  }
  # Define density and survival function for each of the parametric families
  if(family != "gppiece"){
  parlist <- .npar_elife(family = family, par = par)
  valid_par <- try(do.call(what = check_elife_dist,
                           args = parlist), silent = TRUE)
  if(inherits(valid_par, "try-error")){
    return(1e20)
  }
  ldensf <- function(par, dat){
    par$log <- TRUE
    par$x <- dat
    do.call(delife, args = par)
  }
  lsurvf <- function(par, dat, lower.tail = FALSE, log.p = TRUE){
    par$lower.tail <- lower.tail
    par$log.p <- log.p
    par$q <- dat
    do.call(pelife, args = par)
  }

  # Check support constraints
 if(family == "gp"){
    maxdat <- max(time)
    if(par[2] < (-1+1e-8)){
      return(1e20)
    }
    if(par[2] < 0 && (-par[1]/par[2] < maxdat)){
      return(1e20)
    }
  } else if(family == "extgp"){
    maxdat <- max(time)
    # both shape parameters can be negative, but the data
    # must lie inside the support of the distribution
    if(abs(par[3]) <= 1e-8 && par[2] <= 1e-8){
      bounds <- par[1]
    } else if(abs(par[3]) > 1e-8 && par[2] == 0){
      bounds <- c(par[1], ifelse(par[3] < 0, 1+par[3]*maxdat/par[1], Inf))
    } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){
      bounds <- par[1:2]
    } else {
      bounds <- c(par[1:2], 1+par[3]*(exp(par[2]*maxdat/par[1])-1)/par[2])
    }
    if(isTRUE(any(bounds < 0))){
      return(1e20)
    }
  } else if(family == "extweibull"){
    maxdat <- max(time)
    if(par[2] < (-1+1e-8)){
      return(1e20)
    }
    if(par[2] < 0 && (1+par[2]*(maxdat/par[1])^(par[3]) < 0)){
      return(1e20)
    }
  }
  par <- parlist
  } else if(family == "gppiece"){
    thresh <- thresh - thresh[1]
    # Check constraints for the model
    scale <- par[1]
    shape <- par[-1]
    m <- length(shape)
    maxdat <- max(time)
    stopifnot("Threshold and shape parameter should be of the same length." = length(thresh) == m)
    w <- as.numeric(diff(thresh))
    sigma <- scale + c(0, cumsum(shape[-m]*w))
    fail <- !isTRUE(all(sigma > 0,
                       shape >= -1,
                       ifelse(shape[-m] < 0, thresh[-1] <= thresh[-m] - sigma[-m]/shape[-m], TRUE),
                       ifelse(shape[m] < 0, maxdat <= thresh[m] - sigma[m]/shape[m], TRUE)))
    if(fail){
      return(1e20)
    }
    ldensf <- function(par, dat){
      dgppiece(x = dat, scale = par[1], shape = par[-1], thresh = thresh, log = TRUE)
    }
    lsurvf <- function(par, dat, lower.tail = FALSE, log.p = TRUE){
      pgppiece(q = dat, scale = par[1], shape = par[-1], thresh = thresh, lower.tail = lower.tail, log.p = log.p)
    }
  }

  if(!is.null(rtrunc1)){
    rtr1_cont <- lsurvf(dat = rtrunc1, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    rtr2_cont <- 1
  }
  if(!is.null(rtrunc2)){
      rtr2_cont <- lsurvf(dat = rtrunc2, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    rtr2_cont <- 1
  }
  if(!is.null(ltrunc1)){
    ltr1_cont <- lsurvf(dat = ltrunc1, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    ltr1_cont <- 0
  }
  if(!is.null(ltrunc2)){
    ltr2_cont <- lsurvf(dat = ltrunc2, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    ltr2_cont <- 0
  }
  # Likelihoods, depending on the sampling scheme (left- and right-truncated, left-truncated and right-censored, none)
  ll <- sum(weights*ldensf(dat = time, par = par) -
                  log(rtr1_cont - ltr1_cont +
                        ifelse(is.na(rtr2_cont), 0, rtr2_cont) -
                        ifelse(is.na(ltr2_cont), 0, ltr2_cont))
            )
  if (!is.finite(ll)) {
    return(1e20)
  }  else {
    return(-ll)
  }
}

#' Fit excess lifetime models for doubly interval truncated data
#'
#' This function is a wrapper around constrained optimization
#' routines for different models with non-informative
#' censoring and truncation patterns.
#'
#' @note The extended generalized Pareto model is constrained
#' to avoid regions where the likelihood is flat so \eqn{\xi \in [-1, 10]} in
#' the optimization algorithm.
#'
#' The standard errors are obtained via the observed information matrix, calculated
#' using the hessian. In many instances, such as when the shape parameter is zero
#' or negative, the hessian is singular and no estimates are returned.
#'
#' @keywords internal
#' @inheritParams nll_ditrunc_elife
#' @inheritParams nll_elife
#' @param export logical; should data be included in the returned object to produce diagnostic plots? Default to \code{FALSE}.
#' @param start vector of starting values for the optimization routine. If \code{NULL}, the algorithm attempts to find default values and returns a warning with
#' false convergence diagnostic if it cannot.
#' @param restart logical; should multiple starting values be attempted? Default to \code{FALSE}.
#' @return an object of class \code{elife_par}
#' @export
fit_ditrunc_elife <- function(
  time,
  ltrunc1 = NULL,
  rtrunc1 = NULL,
  ltrunc2 = NULL,
  rtrunc2 = NULL,
  thresh = 0,
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
  weights = NULL,
  export = FALSE,
  start = NULL,
  restart = FALSE,
  arguments = NULL, ...
  ){
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = fit_ditrunc_elife, call = call, arguments = arguments)
    return(do.call(fit_ditrunc_elife, args = arguments))
  }
  stopifnot("Argument `restart` should be a logical vector" = is.logical(restart) & length(restart) == 1L)
  if(is.null(weights)){
    weights <- rep(1, length(time))
  } else{
    stopifnot("All weights should be positive" = isTRUE(all(weights >=0 )))
  }
  family <- match.arg(family)
  stopifnot("Threshold must be positive" = thresh >= 0,
            "Only a single threshold is allowed" = length(thresh) == 1L | family == "gppiece"
  )
  if(family == "gpppiece" && length(thresh) == 1L){
    family <- "gp"
  }
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- time > thresh[1]
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    if(!is.null(ltrunc1)){ #both ltrc and ltrt
      ltrunc1 <- pmax(0, ltrunc1[ind] - thresh[1])
    }
    if(!is.null(ltrunc2)){
      ltrunc2 <- pmax(0, ltrunc2[ind] - thresh[1])
    }
    if(!is.null(rtrunc1)){
      rtrunc1 <- pmax(0, rtrunc1[ind] - thresh[1])
    }
    if(!is.null(rtrunc2)){
      rtrunc2 <- pmax(0, rtrunc2[ind] - thresh[1])
    }
    # Sometimes, the thresholding leads to ltrunc1=rtrunc1=0
    bad <- which(rtrunc1 == 0)
    if(length(bad) > 0){
      ltrunc1[bad] <- ltrunc2[bad]
      rtrunc1[bad] <- rtrunc2[bad]
      ltrunc2[bad] <- NA
      rtrunc2[bad] <- NA
    }
  }



  # Keep maximum and sample size
  n <- length(time) #number of exceedances
  maxdat <- max(time)
  stopifnot("Incorrect data input: some entries are missing or finite" = is.finite(maxdat))
  # For gppiece model
  # thresh <- thresh - thresh[1]
  # Perform constrained optimization for the different routines
  if(family == "exp"){
    # closed-form solution for left-truncated and right-censored data
    if(is.null(start)){
        start <- mean(time)
      } else{
        stopifnot("`start should be a scalar." = length(start) == 1L,
                  "`The parameter of the exponential distribution should be positive." = start > 0)
      }
      opt_mle <- optim(par = start,
                       fn = nll_ditrunc_elife,
                       method = "Brent",
                       lower = 1e-8,
                       upper = 3 * maxdat,
                       thresh = 0,
                       time = time,
                       ltrunc1 = ltrunc1,
                       rtrunc1 = rtrunc1,
                       ltrunc2 = ltrunc2,
                       rtrunc2 = rtrunc2,
                       weights = weights,
                       family = "exp",
                       hessian = TRUE
      )
      mle <- opt_mle$par
      vcov <- try(solve(opt_mle$hessian))
      if(is.character(vcov)){
        vcov <- NULL
        se_mle <- NA
      } else{
        se_mle <- try(sqrt(diag(vcov)))
      }
      if(is.character(se_mle)){
        se_mle <- NA
      }
      ll <- -opt_mle$value
      conv <- opt_mle$convergence == 0
    } else {
      if(family == "gppiece"){
      m <- length(thresh)
      hin <- function(par, maxdat = NULL, thresh = 0, ...){
        stopifnot("Argument \"maxdat\" is missing" = !is.null(maxdat),
                  "Vector of threshold should not be a single value" = length(thresh) > 1)
        m <- length(thresh)
        w <- as.numeric(diff(thresh))
        shape <- par[-1]
        sigma <- par[1] + c(0, cumsum(shape[-m]*w))
          c(sigma,
          shape,
          ifelse(shape[-m] < 0, thresh[-m] - sigma[-m]/shape[-m] - thresh[-1], 1e-5),
          ifelse(shape[m] < 0, thresh[m] - sigma[m]/shape[m] - maxdat, 1e-5)
        )
      }
      LB <- c(0, rep(-1, m))
      UB <- c(Inf, rep(2, m))
      ineqLB <- c(rep(0, m), rep(-1, m), rep(0, m)) # Constraints for xi...
      ineqUB <- c(rep(Inf, m), rep(4, m), rep(Inf, m)) # Careful, these are for GPD
      if(is.null(start)){
      st <- try(fit_ditrunc_elife(time = time,
                          thresh = 0,
                          ltrunc1 = ltrunc1,
                          rtrunc1 = rtrunc1,
                          ltrunc2 = ltrunc2,
                          rtrunc2 = rtrunc2,
                          family = "gp",
                          weights = weights))
      if(!is.character(st) && st$convergence){
        start <- c(st$par['scale'], rep(st$par['shape'], m))
      } else{
        start <- c(mean(time, na.rm = TRUE), rep(0.04, m))
      }
     } else{
        stopifnot("Incorrect parameter length." = length(start) == (m + 1L))
        ineq <- hin(start, maxdat = maxdat, thresh = thresh-thresh[1])
        stopifnot("Invalid starting values." = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
      }
      } else {
          ineq_fn <- .ineq_optim_elife(dat = time,
                                       family = family,
                                       maxdat = maxdat,
                                       start = start)
          LB <- ineq_fn$LB
          UB <- ineq_fn$UB
          ineqLB <- ineq_fn$ineqLB
          ineqUB <- ineq_fn$ineqUB
          hin <- ineq_fn$hin
          start <- ineq_fn$start
      }

    stopifnot("Invalid starting values." =
        nll_ditrunc_elife(par = start,
                          family = family,
                          ltrunc1 = ltrunc1,
                          ltrunc2 = ltrunc2,
                          rtrunc1 = rtrunc1,
                          rtrunc2 = rtrunc2,
                          time = time,
                          thresh = thresh - thresh[1],
                          weights = weights) < 1e20)
    opt_mle <- Rsolnp::solnp(pars = start,
                             family = family,
                             ltrunc1 = ltrunc1,
                             ltrunc2 = ltrunc2,
                             rtrunc1 = rtrunc1,
                             rtrunc2 = rtrunc2,
                             time = time,
                             fun = nll_ditrunc_elife,
                             ineqfun = hin,
                             maxdat = maxdat,
                             thresh = thresh - thresh[1],
                             weights = weights,
                             ineqLB = ineqLB,
                             ineqUB = ineqUB,
                             control = list(trace = 0))
    if(isTRUE(opt_mle$convergence != 0 | restart)){
      opt_mle <- Rsolnp::gosolnp(LB = LB,
                                 UB = ifelse(is.finite(UB), UB, 10*maxdat),
                                 family = family,
                                 ltrunc1 = ltrunc1,
                                 ltrunc2 = ltrunc2,
                                 rtrunc1 = rtrunc1,
                                 rtrunc2 = rtrunc2,
                                 time = time,
                                 fun = nll_ditrunc_elife,
                                 ineqfun = hin,
                                 maxdat = maxdat,
                                 thresh = thresh - thresh[1],
                                 weights = weights,
                                 ineqLB = ineqLB,
                                 ineqUB = ineqUB,
                               control = list(trace = 0),
                               n.sim = 200L,
                               n.restarts = 10L)
    }
    mle <- opt_mle$pars
    vcov <- try(solve(opt_mle$hessian[-(seq_along(ineqLB)),-(seq_along(ineqLB))]))
    if(is.character(vcov)){
      vcov <- NULL
      se_mle <- rep(NA, length(start))
    } else{
      se_mle <- try(sqrt(diag(vcov)))
    }
    if(is.character(se_mle)){
      se_mle <- rep(NA, length(start))
    }
    # The function returns -ll, and the value at each
    # iteration (thus keep only the last one)
    ll <- -opt_mle$values[length(opt_mle$values)]
    conv <- opt_mle$convergence == 0
  }
names(mle) <- names(se_mle) <-
  switch(family,
         "exp" = c("scale"),
         "gomp" = c("scale","shape"),
         "gompmake" = c("scale","lambda","shape"),
         "weibull" = c("scale","shape"),
         "extweibull" = c("scale","shape", "xi"),
         "gp" = c("scale","shape"),
         "extgp" = c("scale","beta","gamma"),
         "gppiece" = c("scale", paste0("shape",1:(length(mle)-1L))),
         "perks" = c("rate","shape"),
         "perksmake" = c("rate","lambda","shape"),
         "beard" = c("rate","alpha","beta"),
         "beardmake" = c("rate","lambda","alpha","beta")
  )
    cens_type <- "none"

    if(!isTRUE(all(is.na(rtrunc2), is.na(ltrunc2)))){
      trunc_type <- "doubly interval truncated"
    } else{
      trunc_type <- "interval truncated"
    }
    if(ll == -1e20){
    conv <- FALSE
    warning("Algorithm did not converge: try changing the starting values.")
    }
  if(isTRUE(export)){

    ret <- structure(list(par = mle,
                   std.error = se_mle,
                   loglik = ll,
                   nexc = sum(weights),
                   vcov = vcov,
                   convergence = conv,
                   family = family,
                   thresh = thresh,
                   time = time,
                   weights = weights,
                   ltrunc = cbind(ltrunc1, ltrunc2),
                   rtrunc = cbind(rtrunc1, rtrunc2),
                   cens_type = cens_type,
                   trunc_type = trunc_type),
              class = "elife_par")
    return(ret)
  } else{
    structure(list(par = mle,
                   std.error = se_mle,
                   loglik = ll,
                   nexc = sum(weights),
                   vcov = vcov,
                   convergence = conv,
                   cens_type = cens_type,
                   trunc_type = trunc_type,
                   family = family,
                   thresh = thresh),
              class = "elife_par")
  }
}

#' Likelihood ratio test for doubly interval truncated data
#'
#' This function fits separate models for each distinct
#' value of covariates and computes a likelihood ratio test
#' to test whether there are significant differences between
#' groups.
#'
#' @keywords internal
#' @export
#' @inheritParams nll_elife
#' @param covariate vector of factors, logical or integer whose distinct values are
#' @return a list with elements
#' \itemize{
#' \item \code{stat}: likelihood ratio statistic
#' \item \code{df}: degrees of freedom
#' \item \code{pval}: the p-value obtained from the asymptotic chi-square approximation.
#' }
test_ditrunc_elife <- function(time,
                       covariate,
                       thresh = 0,
                       ltrunc1 = NULL,
                       rtrunc1 = NULL,
                       ltrunc2 = NULL,
                       rtrunc2 = NULL,
                       family = c("exp",
                                  "gp",
                                  "gomp",
                                  "gompmake",
                                  "weibull",
                                  "extgp",
                                  "extweibull",
                                  "perks",
                                  "beard",
                                  "perksmake",
                                  "beardmake"),
                       weights = rep(1, length(time)),
                       arguments = NULL, ...) {
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = test_ditrunc_elife, call = call, arguments = arguments)
    return(do.call(test_ditrunc_elife, args = arguments))
  }
  family <- match.arg(family)
  stopifnot("Covariate must be provided" = !missing(covariate),
            "Object `covariate` should be of the same length as `time`" = length(covariate) == length(time),
            "Provide a single threshold" = !missing(thresh) && length(thresh) == 1L)
  npar <- switch(family,
                 "exp" = 1L,
                 "gp" = 2L,
                 "gomp" = 2L,
                 "gompmake" = 3L,
                 "extgp" = 3L,
                 "weibull" = 2L,
                 "extweibull" = 3L,
                 "perks" = 2L,
                 "beard" = 3L,
                 "perksmake" = 3L,
                 "beardmake" = 4L)
  # Transform to factor
  covariate <- as.factor(covariate)
  nobs_cov <- table(covariate)
  m <- length(nobs_cov)
  stopifnot("There should be more than one group in `covariate`." = m > 1,
            "There are too few observations (less than 5 times the number of parameters) for some modalities of `covariate`." = min(nobs_cov) >= 5*npar)
  # Fit the pooled model
  labels <- names(nobs_cov)
  fit_null <- try(fit_ditrunc_elife(time = time,
                            thresh = thresh,
                            ltrunc1 = ltrunc1,
                            rtrunc1 = rtrunc1,
                            ltrunc2 = ltrunc2,
                            rtrunc2 = rtrunc2,
                            family = family,
                            weights = weights))
  loglik0 <- ifelse(is.character(fit_null), NA, fit_null$loglik)
  fit_alternative <- list()
  loglik1 <- rep(0, m)
  n_levels <- rep(0L, m)
  for(i in 1:m){
    fit_alternative[[i]] <- try(fit_ditrunc_elife(time = time[covariate == labels[i]],
                                          thresh = thresh,
                                          ltrunc1 = ltrunc1[covariate == labels[i]],
                                          rtrunc1 = rtrunc1[covariate == labels[i]],
                                          ltrunc2 = ltrunc2[covariate == labels[i]],
                                          rtrunc2 = rtrunc2[covariate == labels[i]],
                                          family = family,
                                          weights = weights[covariate == labels[i]]))
    loglik1[i] <- ifelse(is.character(fit_alternative[[i]]), NA, fit_alternative[[i]]$loglik)
    n_levels[i] <- fit_alternative[[i]]$nexc
  }
  lrt_stat <- 2*as.numeric((sum(loglik1)-loglik0))
  names(n_levels) <- labels
  p_value <- pchisq(q = lrt_stat, df = (m - 1) * npar, lower.tail = FALSE)
  invisible(structure(list(stat = lrt_stat,
                           df = (m - 1) * npar,
                           pval = p_value,
                           nobs_covariate = n_levels,
                           thresh = thresh,
                           family = family),
                      class = "elife_par_test"))
}
