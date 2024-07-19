#' Likelihood for arbitrary censored and truncated data
#'
#' Computes the log-likelihood for various parametric models suitable for threshold
#' exceedances. If threshold is non-zero, then only right-censored, observed event time and interval censored
#' data whose timing exceeds the thresholds are kept.
#' @param par vector of parameters, in the following order: scale, rate and shape
#' @param thresh vector of thresholds
#' @inheritParams npsurv
#' @param ltrunc lower truncation limit, default to \code{NULL}
#' @param rtrunc upper truncation limit, default to \code{NULL}
#' @param weights weights for observations
#' @param family string; choice of parametric family
#' @param status integer vector giving status of an observation. If \code{NULL} (default), this argument is computed internally based on \code{type}.
#' @param ... additional arguments for optimization, currently ignored.
#' @return log-likelihood values
#' @export
#' @examples
#' data(ewsim, package = "longevity")
#' nll_elife(par = c(5, 0.3),
#'           family = "gp",
#'           arguments = ewsim)
nll_elife <- function(par,
                      time,
                      time2 = NULL,
                      event = NULL,
                      type = c("right",
                               "left",
                               "interval",
                               "interval2"),
                      ltrunc = NULL,
                      rtrunc = NULL,
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
                      weights = NULL,
                      status = NULL,
                      arguments = NULL,
                      ...){
  if(missing(time)){
    time <- NULL
  }
  if(missing(par)){
    par <- NULL
  }
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = "nll_elife", call = call, arguments = arguments)
    return(do.call(nll_elife, args = arguments))
  }
  if(is.null(time)){
    stop("argument \"time\" is missing, with no default")
  }
  if(is.null(par)){
    stop("argument \"par\" is missing, with no default")
  }
  if(is.null(weights)){
    weights <- rep(1, length(time))
  }
  family <- match.arg(family)
  type <- match.arg(type)
  # All parameters must be finite values
  stopifnot(isTRUE(all(is.finite(par))))
  stopifnot("Threshold must be non-negative" = thresh >= 0,
            "Only a single threshold value is allowed" = family == "gppiece" || length(thresh) == 1L
  )
  if(isTRUE(all(is.matrix(ltrunc),
         is.matrix(rtrunc),
         ncol(ltrunc) == ncol(rtrunc),
         ncol(rtrunc) == 2L))){
    # Doubly truncated data
    stopifnot("Censoring is not currently handled for doubly truncated data." = is.null(event) | isTRUE(all(event == 1L)),
              "Argument `time2` not used for doubly truncated data" = is.null(time2)
              )
    nll_ditrunc_elife(
      par = par,
      time = time,
      ltrunc1 = ltrunc[,1],
      rtrunc1 = rtrunc[,1],
      ltrunc2 = ltrunc[,2],
      rtrunc2 = rtrunc[,2],
      family = family,
      thresh = thresh,
      weights = weights)

  } else{
  # Format time entry and check for input
  if(is.null(status)){
    survout <- .check_surv(time = time,
                           time2 = time2,
                           event = event,
                           type = type)
    time <- survout$time
    time2 <- survout$time2
    status <- survout$status
  }
  stopifnot("Status could not be resolved." = isTRUE(all(status %in% 0:3)),
            "Incorrect `ltrunc` argument." = is.null(ltrunc) || (!is.null(ltrunc) && (length(ltrunc) %in% c(1L, length(time)))),
            "Incorrect `rtrunc` argument." = is.null(rtrunc) || (!is.null(rtrunc) && (length(rtrunc) %in% c(1L, length(time)))))
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2, FALSE, ifelse(status == 3, time >= thresh[1], time > thresh[1]))
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    time2 <- time2[ind] - thresh[1]
    status <- status[ind]
    if(!is.null(ltrunc)){
      ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
    }
    if(!is.null(rtrunc)){
      rtrunc <- rtrunc[ind] - thresh[1]
    }
  }
  if(!is.null(ltrunc)){
    if(isTRUE(any(ltrunc > time, ltrunc > time2, na.rm = TRUE))){
      stop("Left-truncation must be lower than observation times.")
    }
  }
  if(!is.null(rtrunc)){
    if(isTRUE(any(rtrunc < time, rtrunc < time2, na.rm = TRUE))){
      stop("Right-truncation must be lower than observation times.")
    }
  }
  if(family != "gppiece"){
  parlist <- .npar_elife(family = family, par = par)
  valid_par <- try(do.call(what = check_elife_dist,
                           args = parlist), silent = TRUE)
  if(inherits(valid_par, "try-error")){
    return(1e20)
  }
  ldensf <- function(par, dat){
    stopifnot(is.list(par))
    par$log <- TRUE
    par$x <- dat
    do.call(delife, args = par)
  }
  lsurvf <- function(par, dat, lower.tail = FALSE, log.p = TRUE){
    stopifnot(is.list(par))
    par$lower.tail <- lower.tail
    par$log.p <- log.p
    par$q <- dat
    do.call(pelife, args = par)
  }

  # Check support constraints
 if(family == "gp"){
    maxdat <- max(ifelse(status == 2L, time2, time))
    if(par[2] < (-1+1e-8)){
      return(1e20)
    }
    if(par[2] < 0 && (-par[1]/par[2] < maxdat)){
      return(1e20)
    }
  } else if(family == "extgp"){
    maxdat <- max(ifelse(status == 2L, time2, time))
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
  # Replace parameters by list for call
  par <- parlist
  } else if(family == "gppiece"){
    stopifnot("Length of parameters for piecewise generalized Pareto \nfamily is the number of threshold values plus one." = length(par) == length(thresh) + 1L)
    thresh <- thresh - thresh[1]
    # Check constraints for the model
    scale <- par[1]
    shape <- par[-1]
    m <- length(shape)
    maxdat <- max(ifelse(status == 2L, time2, time))
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
  if(!is.null(rtrunc)){
    rtr_cont <- lsurvf(dat = rtrunc, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    rtr_cont <- 1
  }
  if(!is.null(ltrunc)){
    ltr_cont <- lsurvf(dat = ltrunc, par = par, lower.tail = TRUE, log.p = FALSE)
  } else{
    ltr_cont <- 0
  }
  # Likelihoods, depending on the sampling scheme (left- and right-truncated, left-truncated and right-censored, none)
  ll <- sum(weights*(ifelse(status == 1L, ldensf(dat = time, par = par),
                            ifelse(status == 0L, log(pmax(0, rtr_cont - lsurvf(dat = time, par = par, lower.tail = TRUE, log.p = FALSE))),
                                   ifelse(status == 2L, log(pmax(0, lsurvf(dat = time2, par = par, lower.tail = TRUE, log.p = FALSE) - ltr_cont)),
                                          log(pmin(rtr_cont, lsurvf(dat = time2, par = par, lower.tail = TRUE, log.p = FALSE)) - pmax(ltr_cont, lsurvf(dat = time, par = par, lower.tail = TRUE, log.p = FALSE)))))) -
                       log(rtr_cont - ltr_cont)))
  if (!is.finite(ll)) {
    return(1e20)
  }  else {
    return(-ll)
  }
  }
}

.npar_elife <- function(par, family, return_npar = FALSE){
  stopifnot(is.logical(return_npar))
  family <- match.arg(family, choices =  c("exp", "gp", "gomp", "gompmake", "weibull", "extgp", "gppiece",
                                           "extweibull", "perks", "beard", "perksmake", "beardmake"))
   # Number of scale, rate and shape parameters
  npars <- switch(family,
                  "exp" = c(1, 0, 0),
                  "gp" = c(1, 0, 1),
                  "weibull" = c(1, 0, 1),
                  "gomp" = c(1, 0, 1),
                  "gompmake" = rep(1, 3),
                  "extgp" = c(1, 0, 2),
                  "gppiece" = c(1, 0, length(par) - 1),
                  "extweibull" = c(1, 0, 2),
                  "perks" = c(0, 1, 1),
                  "beard" = c(0, 1, 2),
                  "perksmake" = c(0, 2, 1),
                  "beardmake" = c(0, 2, 2))
  if(isTRUE(return_npar)){
    return(npars)
  }
  if(family == "gppiece" & length(par) < 2){
    stop("Invalid parameter vector for model \"gppiece\".")
  }
  if(length(par) != sum(npars)){
    stop(paste0("Invalid parameter length for family \"", family, "\"."))
  }
  results <- list(family = family)
  if(npars[1] > 0){
    results$scale <- par[seq_len(npars[1])]
  }
  if(npars[2] > 0){
    results$rate <- par[seq_len(npars[2]) + npars[1]]
  }
  if(npars[3] > 0){
    results$shape <- par[seq_len(npars[3]) + sum(npars[1:2])]
  }
  return(results)
}

#' Fit excess lifetime models by maximum likelihood
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
#' @importFrom Rsolnp solnp
#' @inheritParams nll_elife
#' @param export logical; should data be included in the returned object to produce diagnostic plots? Default to \code{FALSE}.
#' @param start vector of starting values for the optimization routine. If \code{NULL}, the algorithm attempts to find default values and returns a warning with
#' false convergence diagnostic if it cannot.
#' @param restart logical; should multiple starting values be attempted? Default to \code{FALSE}.
#' @param arguments a named list specifying default arguments of the function that are common to all \code{elife} calls
#' @param check logical; if \code{TRUE}, fit all submodels to ensure that simpler models fit worst or as well
#' @param ... additional parameters, currently ignored
#' @return an object of class \code{elife_par}
#' @export
#' @examples
#' data(ewsim, package = "longevity")
#' fit1 <- fit_elife(
#'    arguments = ewsim,
#'    export = TRUE,
#'    family = "exp")
#' fit2 <- fit_elife(
#'    arguments = ewsim,
#'    export = TRUE,
#'    family = "gp")
#' plot(fit1)
#' summary(fit1)
#' anova(fit2, fit1)
fit_elife <- function(time,
                      time2 = NULL,
                      event = NULL,
                      type = c("right","left","interval","interval2"),
                      ltrunc = NULL,
                      rtrunc = NULL,
                      thresh = 0,
                      status = NULL,
                      family = c("exp", "gp", "weibull", "gomp",
                                 "gompmake", "extgp", "gppiece", "extweibull",
                                 "perks", "perksmake", "beard", "beardmake"),
                      weights = NULL,
                      export = FALSE,
                      start = NULL,
                      restart = FALSE,
                      arguments = NULL,
                      check = FALSE,
                      ...
){
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = fit_elife, call = call, arguments = arguments)
    return(do.call(fit_elife, args = arguments))
  }

  stopifnot("Argument `restart` should be a logical vector" = is.logical(restart) & length(restart) == 1L)
  stopifnot("Argument `event` must be NULL, a vector of the same length as time or a scalar." = isTRUE(is.null(event) | length(event) %in% c(1L, length(time))))
  family <- match.arg(family)
  if(!is.null(event) & length(event) == 1L){
    event <- rep(event, length.out = length(time))
  }
  if(isTRUE(all(is.matrix(ltrunc),
                is.matrix(rtrunc),
                ncol(ltrunc) == ncol(rtrunc),
                ncol(rtrunc) == 2L))){
    # Doubly truncated data
    stopifnot("Censoring is not currently handled for doubly truncated data." = is.null(event) | isTRUE(all(event == 1L)),
              "Argument `time2` not used for doubly truncated data" = is.null(time2)
    )
    return(fit_ditrunc_elife(
      time = time,
      ltrunc1 = ltrunc[,1],
      rtrunc1 = rtrunc[,1],
      ltrunc2 = ltrunc[,2],
      rtrunc2 = rtrunc[,2],
      family = family,
      thresh = thresh,
      weights = weights,
      export = export,
      start = start,
      restart = restart))
  } else{
  if(is.null(status)){
    survout <- .check_surv(time = time,
                           time2 = time2,
                           event = event,
                           type = type)
    time <- survout$time
    time2 <- survout$time2
    status <- survout$status
  }
  n <- length(time)
  if(is.null(weights)){
    weights <- rep(1, length(time))
  } else{
    stopifnot("All weights should be positive" = isTRUE(all(weights >=0)))
  }
  type <- match.arg(type)
  stopifnot("Threshold must be positive" = thresh >= 0,
            "Only a single threshold is allowed" = length(thresh) == 1L || family == "gppiece",
            "Only a single sampling scheme in `type` is allowed" = length(type) == 1L
  )
  if(family == "gppiece" && length(thresh) == 1L){
    family <- "gp"
  }
  if(!is.null(ltrunc)){
    stopifnot("`time` and `ltrunc` must be of the same length." = n == length(ltrunc),
              "`ltrunc` must be smaller than `time2`" = isTRUE(all(ltrunc <= time2, na.rm = TRUE)),
              "`ltrunc` must be smaller than `time`" = isTRUE(all(ltrunc <= time, na.rm = TRUE))
    )
  }
  if(!is.null(rtrunc)){
    stopifnot("`time` and `rtrunc` must be of the same length." = n == length(rtrunc),
              "`rtrunc` must be higher than `time2`" = isTRUE(all(rtrunc >= time2, na.rm = TRUE)),
              "`rtrunc` must be higher than `time`" = isTRUE(all(rtrunc >= time, na.rm = TRUE))
    )
  }
  if(!is.null(ltrunc) && !is.null(rtrunc)){
    stopifnot("`ltrunc` must be larger than `ltrunc" = isTRUE(all(rtrunc > ltrunc)))
  }
  if(thresh[1] > 0){
    thresh0 <- thresh
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2, FALSE, ifelse(status == 3, time >= thresh[1], time > thresh[1]))
    weights <- weights[ind] # in this order
    time <- pmax(0, time[ind] - thresh[1])
    time2 <- pmax(0, time2[ind] - thresh[1])
    status <- status[ind]
    if(!is.null(event)){
      event <- event[ind]
    }
    if(!is.null(ltrunc)){ #both ltrc and ltrt
      ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
    }
    if(!is.null(rtrunc)){
      rtrunc <- rtrunc[ind] - thresh[1]
    }
  }
  # Keep maximum and sample size
  n <- length(time) #number of exceedances
  if(n == 0L){
    stop("No observation lies above the threshold.")
  }
  dat <- ifelse(status == 2L, time2, time)
  maxdat <- max(dat)
  stopifnot("Incorrect data input: some entries are missing or infinite" = is.finite(maxdat))
  # For gppiece model
  # thresh <- thresh - thresh[1]
  # Perform constrained optimization for the different routines
  if(family == "exp"){
    # closed-form solution for left-truncated and right-censored data
    if(is.null(ltrunc) && is.null(rtrunc) && isTRUE(all(status == 1L))){
      mle <- weighted.mean(x = time, w = weights, na.rm = TRUE)
      se_mle <- mle / sqrt(sum(weights))
      ll <- sum(weights*dexp(time, rate = 1/mle, log = TRUE))
      conv <- TRUE
    } else if(is.null(rtrunc) & isTRUE(all(status %in% c(0L, 1L))) & isTRUE(sum(status) > 0)){
      if(is.null(ltrunc)){
        ltrunc <- 0
      }
      exp_mle_ltrc <- function(dat, ltrunc, status, weights){
        sum(weights*(dat - ltrunc)) / sum(weights[status == 1L])
      }
      mle <- exp_mle_ltrc(dat = time,
                          ltrunc = ltrunc,
                          status = status,
                          weights = weights)
      se_mle <- mle/sqrt(sum(weights[status == 1L]))
      ll <- - nll_elife(par = mle,
                        time = time,
                        time2 = time2,
                        event = event,
                        status = status,
                        ltrunc = ltrunc,
                        thresh = 0,
                        family = "exp",
                        weights = weights)
      conv <- TRUE
    } else{
      if(is.null(start)){
        start <- mean(dat)
      } else{
        stopifnot("`start should be a scalar." = length(start) == 1L,
                  "`The parameter of the exponential distribution should be positive." = start > 0)
      }
      opt_mle <- optim(par = start,
                       fn = nll_elife,
                       method = "Brent",
                       lower = 1e-8,
                       upper = 3 * maxdat,
                       type = type,
                       thresh = 0,
                       time = time,
                       time2 = time2,
                       status = status,
                       event = event,
                       ltrunc = ltrunc,
                       rtrunc = rtrunc,
                       weights = weights,
                       family = "exp",
                       hessian = TRUE
      )
      mle <- opt_mle$par
      vcov <- try(solve(opt_mle$hessian), silent = TRUE)
      if(inherits(vcov, "try-error")){
        vcov <- NULL
        se_mle <- NA
      } else{
      se_mle <- try(sqrt(diag(vcov)), silent = TRUE)
      if(is.character(se_mle)){
        se_mle <- rep(NA, length(mle))
      }
      }
      ll <- -opt_mle$value
      conv <- opt_mle$convergence == 0
    }
  } else{
  if(family == "gppiece"){
    m <- length(thresh)
    hin <- function(par, maxdat = NULL, thresh = 0, ...){
      stopifnot("Argument \"maxdat\" is missing" = !is.null(maxdat),
                "Vector of threshold should not be a single value" = length(thresh) > 1L)

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
      st <- try(fit_elife(time = time,
                          status = status,
                          time2 = time2,
                          thresh = 0,
                          ltrunc = ltrunc,
                          rtrunc = rtrunc,
                          type = type,
                          family = "gp",
                          weights= weights))
      if(!is.character(st) && st$convergence){
        start <- c(st$par['scale'], rep(st$par['shape'], m))
      } else{
        start <- c(mean(dat, na.rm = TRUE), rep(0.04, m))
      }
    } else{
      stopifnot("Incorrect parameter length." = length(start) == (m + 1L))
      ineq <- hin(start, maxdat = maxdat, thresh = thresh-thresh[1])
      stopifnot("Invalid starting values." = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  } else {
    ineq_fn <- .ineq_optim_elife(dat = dat,
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
                nll_elife(par = start,
                          family = family,
                          ltrunc = ltrunc,
                          rtrunc = rtrunc,
                          time = time,
                          time2 = time2,
                          status = status,
                          thresh = thresh - thresh[1],
                          type = type,
                          weights = weights) < 1e20)
    opt_mle <- Rsolnp::solnp(pars = start,
                             family = family,
                             ltrunc = ltrunc,
                             rtrunc = rtrunc,
                             time = time,
                             time2 = time2,
                             status = status,
                             fun = nll_elife,
                             ineqfun = hin,
                             maxdat = maxdat,
                             thresh = thresh - thresh[1],
                             type = type,
                             weights = weights,
                             ineqLB = ineqLB,
                             ineqUB = ineqUB,
                             control = list(trace = 0))
    if(opt_mle$convergence != 0 & restart){
      opt_mle <- Rsolnp::gosolnp(LB = LB,
                                 UB = ifelse(is.finite(UB),
                                             UB, 10*maxdat),
                                 family = family,
                                 ltrunc = ltrunc,
                                 rtrunc = rtrunc,
                                 time = time,
                                 time2 = time2,
                                 status = status,
                                 fun = nll_elife,
                                 ineqfun = hin,
                                 maxdat = maxdat,
                                 thresh = thresh - thresh[1],
                                 type = type,
                                 weights = weights,
                                 ineqLB = ineqLB,
                                 ineqUB = ineqUB,
                                 control = list(trace = 0),
                                 n.sim = 200L,
                                 n.restarts = 10L)
    }
    mle <- opt_mle$pars
    vcov <- try(solve(opt_mle$hessian[-(seq_along(ineqLB)),-(seq_along(ineqLB))]))
    if(inherits(vcov, "try-error")){
      vcov <- NULL
      se_mle <- rep(NA, length(mle))
    } else{
      se_mle <- try(sqrt(diag(vcov)), silent = TRUE)
    if(is.character(se_mle)){
      se_mle <- rep(NA, length(mle))
    }
    }
    if(is.character(se_mle)){
      se_mle <- rep(NA, length(mle))
    }
    # The function returns -ll, and the value at each
    # iteration (thus keep only the last one)
    ll <- -opt_mle$values[length(opt_mle$values)]
    conv <- opt_mle$convergence == 0
  }
  }
  names(mle) <- names(se_mle) <-
    switch(family,
          "exp" = c("scale"),
          "gomp" = c("scale","shape"),
          "gompmake" = c("scale","lambda","shape"),
          "weibull" = c("scale","shape"),
          "extweibull" = c("scale","shape", "xi"),
          "gp" = c("scale","shape"),
          "extgp" = c("scale","beta","xi"),
          "gppiece" = c("scale", paste0("shape",1:(length(mle)-1L))),
          "perks" = c("rate","shape"),
          "perksmake" = c("rate","lambda","shape"),
          "beard" = c("rate","alpha","beta"),
          "beardmake" = c("rate","lambda","alpha","beta")
  )
  if(!is.null(vcov) & is.matrix(vcov)){
    # Name columns of vcov for confint.default routine
    colnames(vcov) <- rownames(vcov) <- names(mle)
  }
  if(isTRUE(all(status == 1L, na.rm = TRUE))){
    cens_type <- "none"
  } else if(type == "right" || isTRUE(all(status %in% 0:1))){
    cens_type <- "right censored"
  } else if(type == "left" || isTRUE(all(status %in% 1:2))){
    cens_type <- "left censored"
  } else{
    cens_type <- "interval censored"
  }
  if(is.null(ltrunc) && is.null(rtrunc)){
    trunc_type <- "none"
  } else if(!is.null(ltrunc) && is.null(rtrunc)){
    trunc_type <- "left truncated"
  } else if(is.null(ltrunc) && !is.null(rtrunc)){
    trunc_type <- "right truncated"
  } else{
    trunc_type <- "interval truncated"
  }
  if(ll == -1e20){
    conv <- FALSE
    warning("Algorithm did not converge: try changing the starting values.")
  }
  # Check that model is "optimal"
  if(isTRUE(check[1]) & family != "exp"){
    sub <- "exp"
    if(family %in% c("extgp","gppiece","extweibull")){
      sub <- c(sub, "gp")
    }
    if(family == "extweibull"){
     sub <- c(sub, "weibull")
    }
    if(family %in% c("extgp","gompmake","beardmake","beard")){
     sub <- c(sub, "gomp")
    }
    if(family == "beardmake"){
     sub <- c(sub, "gompmake","beard","perksmake")
    }
    if(family %in% c("beard","perksmake")){
     sub <- c(sub, "perks")
    }
  # Call all submodels and check fit
  fits <- list()
  for(fam in seq_along(sub)){
    fits[[fam]] <- fit_elife(
      time = time,
      time2 = time2,
      event = event,
      family = sub[fam],
      type = type,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      thresh = thresh,
      status = status,
      weights = weights,
      export = FALSE,
      restart = restart,
      check = FALSE)
  }
  # Extract log likelihood
  devs <- sapply(fits, deviance)
  best <- which.min(devs)
  # Check if submodel fits better
  if(isTRUE(min(devs) < -2*ll)){
    # Submodel is better
  bestsubfam <- fits[[best]]$family # extract family
  # Match parameters
  pars <- .matchpars(family0 = bestsubfam, family1 = family, coef = fits[[best]]$par)
  # Call again function with optimal submodel as starting value
  newopt <- fit_elife(
    time = time,
    time2 = time2,
    event = event,
    family = family,
    type = type,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    thresh = thresh,
    status = status,
    weights = weights,
    export = FALSE,
    start = pars$par,
    restart = FALSE,
    check = FALSE)
  if(isTRUE(deviance(newopt) < devs[best])){
   return(newopt)
  }
  # The submodel fits best
  mle <- pars$par
  # Compute negative log likelihood and hessian manually
  llfn <- function(parv){
    - nll_elife(time = time,
                time2 = time2,
                event = event,
                family = family,
                type = type,
                ltrunc = ltrunc,
                rtrunc = rtrunc,
                thresh = thresh,
                status = status,
                weights = weights,
                par = parv)
  }
  mle <- pars$par
  ll <- llfn(parv = mle)
  if(isTRUE(pars$regular)){
    if(!(family %in% c("gp","extgp")) | mle[length(mle)] > -0.5){
    hessian <- numDeriv::hessian(func = llfn, x = mle)
    vcov <- try(solve(hessian))
    }
  } else{
    vcov <- NULL
    se_mle <- rep(NA, length(mle))
  }
   if(inherits(vcov, "try-error")){
     vcov <- NULL
     se_mle <- rep(NA, length(mle))
   } else{
     se_mle <- try(sqrt(diag(vcov)), silent = TRUE)
   }
   if(is.character(se_mle)){
     se_mle <- rep(NA, length(mle))
   }
  if(length(se_mle) == 0L){
    se_mle <- rep(NA, length(mle))
  }
   # The function returns -ll, and the value at each
   # iteration (thus keep only the last one)
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
           "extgp" = c("scale","beta","xi"),
           "gppiece" = c("scale", paste0("shape",1:(length(mle)-1L))),
           "perks" = c("rate","shape"),
           "perksmake" = c("rate","lambda","shape"),
           "beard" = c("rate","alpha","beta"),
           "beardmake" = c("rate","lambda","alpha","beta")
    )
  if(!is.null(vcov) & (is.matrix(vcov) & isTRUE(nrow(vcov) == length(mle)))){
    # Name columns of vcov for confint.default routine
    colnames(vcov) <- rownames(vcov) <- names(mle)
  }

  }
  if(family %in% c("extgp","gp")){
    if(mle[length(mle)] > 1.999){
      conv <- FALSE
    }
  }
  if(isTRUE(export)){
    ret <- structure(list(par = mle,
                          std.error = se_mle,
                          loglik = ll,
                          nexc = sum(weights),
                          vcov = vcov,
                          convergence = conv,
                          type = type,
                          family = family,
                          thresh = thresh,
                          time = time,
                          time2 = time2,
                          event = event,
                          status = status,
                          weights = weights,
                          ltrunc = ltrunc,
                          rtrunc = rtrunc,
                          cens_type = cens_type,
                          trunc_type = trunc_type),
                     class = "elife_par")
     if(!type  %in% c("interval","interval2")){
       ret$time2 <- NULL
     }
    return(ret)
  } else{
    structure(list(par = mle,
                   std.error = se_mle,
                   loglik = ll,
                   nexc = sum(weights),
                   vcov = vcov,
                   convergence = conv,
                   type = type,
                   cens_type = cens_type,
                   trunc_type = trunc_type,
                   family = family,
                   thresh = thresh),
              class = "elife_par")
  }
}

.matchpars <- function(coef, family0, family1, thresh){
  family0 <- match.arg(arg = family0,
                       choices = c("exp", "gp", "weibull", "gomp",
                                   "gompmake", "extgp", "gppiece", "extweibull",
                                   "perks", "perksmake", "beard", "beardmake"),
                       several.ok = FALSE)
  family1 <- match.arg(arg = family1,
                       choices = c("exp", "gp", "weibull", "gomp",
                                   "gompmake", "extgp", "gppiece", "extweibull",
                                   "perks", "perksmake", "beard", "beardmake"),
                       several.ok = FALSE)
  stopifnot(is.numeric(coef),
            length(coef) == sum(.npar_elife(family = family0, return_npar = TRUE)))
  nmods <- rbind(
    c("exp", "weibull", "regular"),        # 1
    c("exp", "gp", "regular"),             # 2
    c("exp", "gppiece", "regular"),        # 3
    c("exp", "gomp", "boundary"),          # 4
    c("exp", "extgp", "boundary"),         # 5
    c("exp", "gompmake", "invalid"),       # 6
    c("gp", "extgp", "boundary"),          # 7
    c("gomp", "extgp", "regular"),         # 8
    c("gp", "gppiece", "regular"),         # 9
    c("gomp", "gompmake", "boundary"),     # 10
    c("weibull","extweibull","regular"),   # 11
    c("gp","extweibull","regular"),        # 12
    c("exp","extweibull","regular"),       # 13
    c("exp", "perks","boundary"),          # 14
    c("exp", "beard","boundary"),          # 15
    c("gomp", "beard","boundary"),         # 16
    c("perks", "beard","regular"),         # 17
    c("perks", "perksmake","boundary"),    # 18
    c("beard", "beardmake","boundary"),    # 19
    c("perks", "beardmake","boundary"),    # 20
    c("perksmake", "beardmake","regular"), # 21
    c("gompmake", "beardmake","boundary"), # 22
    c("gomp", "beardmake","boundary2"),    # 23
    c("exp", "beardmake","invalid"),       # 24
    c("exp", "perksmake","invalid")        # 25
  )
  ind <- intersect(which(family0 == nmods[,1]), which(family1 == nmods[,2]))
  if(length(ind) == 0){
    warning("Trying to constrain two non-nested models.")
    return(NULL)
  }
  isRegular <- nmods[ind,3] == "regular"
  pars <- switch(ind,
                 c(coef, 1), # 1
                 c(coef, 0), # 2
                 c(coef, rep(0, length(thresh))), # 3
                 c(coef, 0), # 4
                 c(coef, rep(0, 2)), # 5
                 c(coef, 0, 0), # 6
                 c(coef[1], 0, coef[2]), # 7
                 c(coef, 0), # 8
                 c(coef, rep(0, length(thresh) - 1L)), # 9
                 c(coef[1], 0, coef[2]), # 10
                 c(coef, 0), # 11
                 c(coef[1], 1, coef[2]), # 12
                 c(coef, 1, 0), # 13
                 c(0, 1/coef[1]), # 14
                 c(0, 1/coef[1], 0), # 15
                 c(coef[2]/coef[1],1/coef[1], 0), # 16
                 c(coef, 1), # 17
                 c(coef[1], 0, coef[2]), # 18
                 c(coef[1], 0, coef[2:3]), # 19
                 c(coef[1], 0, coef[2], 1), # 20
                 c(coef, 1), # 21
                 c(coef[3]/coef[1], coef[2], 1/coef[1], 0), # 22
                 c(coef[2]/coef[1], 0, 1/coef[1], 0), # 23
                 c(1/coef, 0, 0, 0), # 24
                 c(1/coef, 0, 0) # 25
)
 list(par = pars, regular = isRegular)
}


.ineq_optim_elife <- function(dat, family, maxdat = max(dat), start = NULL, thresh){
  if(family == "gp"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){
      stopifnot("Argument \"maxdat\" is missing, with no default value." = !is.null(maxdat))
      # scale > 0, xi > -1, xdat < -xi/sigma if xi < 0
      c(par[1], par[2], ifelse(par[2] < 0, thresh - par[1]/par[2] - maxdat, 1e-5))
    }
    ineqLB <- c(0, -1, 0)
    ineqUB <- c(Inf, 10, Inf)
    LB <- c(0, -1)
    UB <- c(Inf, 2)
    # hardcode upper bound for shape parameter, to prevent optim from giving nonsensical output
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE), 0.1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 2L)
      ineq <- hin(start, maxdat = maxdat)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  } else if(family == "weibull"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){par }
    ineqLB <- LB <- c(0,0)
    ineqUB <- UB <- rep(Inf, 2)
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE), 1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 2L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  } else if(family == "gompmake"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){par }
    ineqLB <- LB <- c(0, 0, 0)
    ineqUB <- UB <- rep(Inf, 3)
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE),0.1, 0.5)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 3L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }

    # If shape1=0, then exponential model (but only the sum "1/par[1]+par[3]" is identifiable)
  } else if(family == "gomp"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par[1:2] }
    ineqLB <- LB <- rep(0, 2)
    ineqUB <- UB <- rep(Inf, 2)
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE), 0.1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 2L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  }  else if(family == "extweibull"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par[c(1,3)] }
    ineqLB <- LB <- rep(0, 2)
    ineqUB <- UB <- rep(Inf, 2)
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE), 0.1, 0.5)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 3L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  }  else if(family == "perks"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par }
    ineqLB <- LB <- rep(0, 2)
    ineqUB <- UB <- rep(Inf, 2)
    if(is.null(start)){
      start <- c(1/mean(dat, na.rm = TRUE), 1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 2L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  }  else if(family == "perksmake"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par }
    ineqLB <- LB <- rep(0, 3)
    ineqUB <- UB <- rep(Inf, 3)
    if(is.null(start)){
      start <- c(1/mean(dat, na.rm = TRUE), 0.1, 1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 3L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  }  else if(family == "beard"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par }
    ineqLB <- LB <- rep(0, 3)
    ineqUB <- UB <- rep(Inf, 3)
    if(is.null(start)){
      start <- c(1/mean(dat, na.rm = TRUE), 1, 0.5)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 3L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  }  else if(family == "beardmake"){
    hin <- function(par, maxdat = NULL, thresh = 0, ...){ par }
    ineqLB <- LB <- rep(0, 4)
    ineqUB <- UB <- rep(Inf, 4)
    if(is.null(start)){
      start <- c(1/mean(dat, na.rm = TRUE), 0.1, 1, 0.5)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 4L)
      ineq <- hin(start)
      stopifnot("Invalid starting values" = isTRUE(all(ineq > ineqLB, ineq < ineqUB)))
    }
  } else if(family == "extgp"){
    #parameters are (1) scale > 0, (2) beta >= 0 (3) gamma
    hin <- function(par, maxdat = NULL, thresh = 0, ...){
      stopifnot("Argument \"maxdat\" is missing, with no default value." = !is.null(maxdat))
      c(par,
        ifelse(par[3] < 0, 1 - par[2]/par[3], 1e-5),
        ifelse(par[3] < 0 & par[2] > 0, thresh + par[1]/par[2]*log(1-par[2]/par[3]) - maxdat, 1e-5),
        ifelse(par[3] < 0 & par[2] == 0, thresh + -par[1]/par[3] - maxdat, 1e-5)
      )
    }
    ineqLB <- c(rep(0, 2), -1, rep(0, 3))
    ineqUB <- c(rep(Inf, 2), 2, rep(Inf, 3))
    LB <- c(0,0,-1)
    UB <- c(rep(Inf, 2), 2)
    if(is.null(start)){
      start <- c(mean(dat, na.rm = TRUE), 0.1, 0.1)
    } else{
      stopifnot("Incorrect parameter length." = length(start) == 3L)
      ineq <- hin(start, maxdat = maxdat)
      stopifnot("Invalid starting values" = isTRUE(all(ineq >= ineqLB, ineq <= ineqUB)))
    }
    #TODO try also fitting the GP/EXP/Gompertz and see which is best?
  } else{
    stop("Family not implemented")
  }
  return(list(LB = LB,
              UB = UB,
              ineqLB = ineqLB,
              ineqUB = ineqUB,
              hin = hin,
              start = start))
}


#' @export
print.elife_par <-
  function(x,
           digits = min(3, getOption("digits")),
           na.print = "", ...){

    cat("Model:", switch(x$family,
                         exp = "exponential",
                         gomp = "Gompertz",
                         gompmake = "Gompertz-Makeham",
                         weibull = "Weibull",
                         extgp = "extended generalized Pareto",
                         gp = "generalized Pareto",
                         gppiece = "piecewise generalized Pareto",
                         extweibull = "extended Weibull",
                         perks = "Perks",
                         perksmake = "Perks-Makeham",
                         beard = "Beard",
                         beardmake = "Beard-Makeham"),
        "distribution.", "\n")
    if(x$cens_type != "none" || x$trunc_type != "none"){
      cat("Sampling: ",
          ifelse(x$cens_type == "none", "", x$cens_type),
          ifelse(x$cens_type != "none" && x$trunc_type != "none", ", ",""),
          ifelse(x$trunc_type == "none", "", x$trunc_type),
          "\n", sep = "")
    }
    cat("Log-likelihood:", round(x$loglik, digits), "\n")

    cat("\nThreshold:", round(x$thresh, digits), "\n")
    cat("Number of exceedances:", x$nexc, "\n")
    cat("\nEstimates\n")
    print.default(format(x$par, digits = digits), print.gap = 2, quote = FALSE, ...)
    if (!is.null(x$std.error)) {
      if(!isTRUE(all(is.na(x$std.error)))){
      cat("\nStandard Errors\n")
      print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
      }
    }
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    invisible(x)
  }


#' @export
summary.elife_par <- function(object, ...){
  print(object, ...)
}

#' @importFrom stats logLik
#' @export
logLik.elife_par <- function(object, ...) {
  val <- object$loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

#' @importFrom stats deviance
#' @export
deviance.elife_par <- function(object, ...) {
  val <- as.numeric(-2*object$loglik)
  return(val)
}


#' @importFrom stats AIC
#' @export
AIC.elife_par <- function(object, ..., k = 2) {
  val <- as.numeric(-2*object$loglik) + k * length(coef(object))
  return(val)
}

#' @importFrom stats BIC
#' @export
BIC.elife_par <- function(object, ..., k = 2) {
  AIC(object, ..., k = log(nobs(object)))
}

#' @importFrom stats nobs
#' @export
nobs.elife_par <- function(object, ...) {
  return(object$nexc)
}


#' @importFrom stats coef
#' @export
coef.elife_par <- function(object, ...) {
  return(object$par)
}

#' @importFrom stats vcov
#' @export
vcov.elife_par <- function(object, ...) {
  return(object$vcov)
}
