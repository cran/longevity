#' Likelihood ratio test for covariates
#'
#' This function fits separate models for each distinct
#' value of the factor \code{covariate} and computes a likelihood ratio test
#' to test whether there are significant differences between
#' groups.
#'
#' @export
#' @inheritParams nll_elife
#' @param covariate vector of factors, logical or integer whose distinct values define groups
#' @return a list with elements
#' \itemize{
#' \item \code{stat} likelihood ratio statistic
#' \item \code{df} degrees of freedom
#' \item \code{pval} the p-value obtained from the asymptotic chi-square approximation.
#' }
#' @examples
#' test <- with(subset(dutch, ndays > 39082),
#'  test_elife(
#'  time = ndays,
#'  thresh = 39082L,
#'  covariate = gender,
#'  ltrunc = ltrunc,
#'  rtrunc = rtrunc,
#'  family = "exp"))
#'  test
test_elife <- function(
  time,
  time2 = NULL,
  event = NULL,
  covariate,
  thresh = 0,
  ltrunc = NULL,
  rtrunc = NULL,
  type = c("right", "left", "interval", "interval2"),
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
  weights = rep(1, length(time)),
  arguments = NULL,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = "test_elife",
      call = call,
      arguments = arguments
    )
    return(do.call(test_elife, args = arguments))
  }

  family <- match.arg(family)
  type <- match.arg(type)
  stopifnot(
    "Covariate must be provided" = !missing(covariate),
    "Object `covariate` should be of the same length as `time`" = length(
      covariate
    ) ==
      length(time),
    "Provide a single threshold" = length(thresh) == 1L
  )
  npar <- switch(
    family,
    "exp" = 1L,
    "gp" = 2L,
    "gomp" = 2L,
    "gompmake" = 3L,
    "extgp" = 3L,
    "weibull" = 2L,
    "extweibull" = 3L,
    "perks" = 2L,
    "perksmake" = 3L,
    "beard" = 3L,
    "beardmake" = 4L
  )
  if (
    isTRUE(all(
      is.matrix(ltrunc),
      is.matrix(rtrunc),
      ncol(ltrunc) == ncol(rtrunc),
      ncol(rtrunc) == 2L
    ))
  ) {
    # Doubly truncated data
    stopifnot(
      "Censoring is not currently handled for doubly truncated data." = is.null(
        event
      ) |
        isTRUE(all(event == 1L)),
      "Argument `time2` not used for doubly truncated data" = is.null(time2)
    )
    return(
      test_ditrunc_elife(
        time = time,
        covariate = covariate,
        thresh = thresh,
        ltrunc1 = ltrunc[, 1],
        ltrunc2 = ltrunc[, 2],
        rtrunc1 = rtrunc[, 1],
        rtrunc2 = rtrunc[, 2],
        family = family,
        weights = weights
      )
    )
  }
  survout <- .check_surv(time = time, time2 = time2, event = event, type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  # Transform to factor
  if (any(is.na(covariate))) {
    stop("Covariate vector should not include missing values.")
  }
  # Cast to factor, remove unused levels
  covariate <- factor(covariate, exclude = NULL)
  # Count occurences
  nobs_cov <- table(covariate)
  m <- length(nobs_cov)
  stopifnot(
    "There should be more than one group in `covariate`." = m > 1,
    "There are too few observations (less than 5 times the number of parameters) for some modalities of `covariate`." = min(
      nobs_cov
    ) >=
      5 * npar
  )
  # Fit the pooled model
  labels <- names(nobs_cov)
  fit_null <- try(fit_elife(
    time = time,
    time2 = time2,
    event = event,
    status = status,
    thresh = thresh,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    type = type,
    family = family,
    weights = weights
  ))
  loglik0 <- ifelse(is.character(fit_null), NA, fit_null$loglik)
  fit_alternative <- list()
  loglik1 <- rep(0, m)
  n_levels <- rep(0L, m)
  for (i in seq_len(m)) {
    fit_alternative[[i]] <-
      try(fit_elife(
        time = time[covariate == labels[i]],
        time2 = time2[covariate == labels[i]],
        status = status[covariate == labels[i]],
        thresh = thresh,
        ltrunc = ltrunc[covariate == labels[i]],
        rtrunc = rtrunc[covariate == labels[i]],
        type = type,
        family = family,
        weights = weights[covariate == labels[i]]
      ))
    loglik1[i] <- ifelse(
      is.character(fit_alternative[[i]]),
      NA,
      fit_alternative[[i]]$loglik
    )
    n_levels[i] <- fit_alternative[[i]]$nexc
  }
  lrt_stat <- 2 * as.numeric((sum(loglik1) - loglik0))
  names(n_levels) <- labels
  p_value <- pchisq(q = lrt_stat, df = (m - 1) * npar, lower.tail = FALSE)
  out <- structure(
    list(
      stat = lrt_stat,
      df = (m - 1) * npar,
      pval = p_value,
      nobs_covariate = n_levels,
      thresh = thresh,
      family = family
    ),
    class = "elife_par_test"
  )
  out
}
#' @export
print.elife_par_test <- function(
  x,
  digits = min(3, getOption("digits")),
  na.print = "",
  ...
) {
  cat(
    "Model:",
    switch(
      x$family,
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
      beardmake = "Beard-Makeham"
    ),
    "distribution.",
    "\n"
  )
  cat("Threshold:", round(x$thresh, digits), "\n")
  cat("Number of exceedances per covariate level:\n")
  print.default(x$nobs_covariate)
  cat("\nLikelihood ratio statistic:", format(x$stat, digits = digits))
  cat(paste0("\nNull distribution: chi-square (", x$df, ")\n"))
  cat("Asymptotic p-value: ", format(x$pval, digits = digits), "\n")
  invisible(x)
}

#' @importFrom stats anova
#' @export
anova.elife_par <- function(object, object2, ..., test = "Chisq") {
  if (any(missing(object), missing(object2))) {
    stop("Two models must be specified.")
  }
  #test <- match.arg(test)
  test <- "Chisq" # only option implemented so far
  model1 <- deparse(substitute(object))
  model2 <- deparse(substitute(object2))
  models <- c(model1, model2)
  narg <- 2L
  for (i in 1:narg) {
    if (!inherits(get(models[i], envir = parent.frame()), "elife_par")) {
      stop("Invalid input: use only with objects of class 'elife_par'.")
    }
  }
  p <- length(models)
  npar <- nobs <- integer(p)
  dev <- thresh <- numeric(p)
  conv <- rep(FALSE, 2L)
  censt <- trunct <- family <- character(p)

  for (i in seq_len(narg)) {
    elifemod <- get(models[i], envir = parent.frame())
    dev[i] <- deviance(elifemod)
    npar[i] <- length(elifemod$par)
    thresh[i] <- elifemod$thresh[1]
    nobs[i] <- elifemod$nexc
    conv[i] <- elifemod$convergence
    family[i] <- elifemod$family
    censt[i] <- elifemod$cens_type
    trunct[i] <- elifemod$trunc_type
  }
  if (npar[1] < npar[2]) {
    dev <- dev[2:1]
    npar <- npar[2:1]
    family <- family[2:1]
  }
  if (!isTRUE(all.equal(thresh[1], thresh[2]))) {
    stop("Invalid arguments: the thresholds should be the same.")
  }
  if (!isTRUE(all(conv))) {
    stop("At least one of the optimization failed to converge.")
  }
  if (!isTRUE(all.equal(nobs[1], nobs[2]))) {
    stop("Invalid arguments: the observations should be the same.")
  }
  if (family[1] == family[2]) {
    stop("Models should be of different families.")
  }
  if (!isTRUE(all(length(unique(censt)) == 1L, length(unique(trunct)) == 1L))) {
    stop("Models have different truncation or censoring schemes.")
  }
  # Cases considered
  nmods <- rbind(
    c("exp", "weibull", "regular"), # chisq(1)
    c("exp", "gp", "regular"), # chisq(1)
    c("exp", "gppiece", "regular"), # chisq(K)
    c("exp", "gomp", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("exp", "extgp", "boundary"), # 0.5 chisq(1) + 0.5 chisq(2)
    c("exp", "gompmake", "invalid"),
    c("gp", "extgp", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("gomp", "extgp", "regular"), # chisq(1)
    c("gp", "gppiece", "regular"), # chisq(K-1)
    c("gomp", "gompmake", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("weibull", "extweibull", "regular"), # chisq(1)
    c("gp", "extweibull", "regular"), # chisq(1)
    c("exp", "extweibull", "regular"), # chisq(2)
    c("exp", "perks", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("exp", "beard", "boundary"), # 0.5 chisq(1) + 0.5 chisq(2)
    c("gomp", "beard", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("perks", "beard", "regular"), # chisq(1)
    c("perks", "perksmake", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("beard", "beardmake", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)u
    c("perks", "beardmake", "boundary"), # 0.5 chisq(0) + 0.5 chisq(1)
    c("perksmake", "beardmake", "regular"),
    c("gompmake", "beardmake", "boundary"),
    c("gomp", "beardmake", "boundary2"), # 0.25 chisq(0) + 0.5 chisq(1) + 0.25 chisq(2)
    c("exp", "beardmake", "invalid"),
    c("exp", "perksmake", "invalid")
  )
  match_family <- which(apply(nmods[, 1:2], 1, function(fam) {
    isTRUE(all(family %in% fam))
  }))
  stopifnot("Invalid input: models are not nested" = length(match_family) == 1L)
  if (nmods[match_family, 3] == "invalid") {
    stop(
      "Models are nested, but parameters are not identifiable.\nThe information matrix is singular."
    )
  }
  df <- -diff(npar)
  dvdiff <- diff(dev)
  if (dvdiff < 0 && dvdiff > -1e-4) {
    # Numerical tolerance for zero
    dvdiff <- 0
  }
  if (dvdiff < 0) {
    stop(
      "The alternative model has a lower likelihood value than the null model, indicating convergence problems."
    )
  }
  # if(test == "Chisq"){
  if (nmods[match_family, 3] == "regular") {
    #regular model
    pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
  } else if (nmods[match_family, 3] == "boundary") {
    pval <- 0.5 *
      pchisq(dvdiff, df = df, lower.tail = FALSE) +
      0.5 * pchisq(dvdiff, df = df - 1L, lower.tail = FALSE)
  } else if (nmods[match_family, 3] == "boundary2") {
    pval <- 0.25 *
      pchisq(dvdiff, df = df, lower.tail = FALSE) +
      0.5 * pchisq(dvdiff, df = df - 1L, lower.tail = FALSE) +
      0.25 * pchisq(dvdiff, df = df - 2L, lower.tail = FALSE)
  }
  table <- data.frame(npar, dev, c(NA, df), c(NA, dvdiff), c(NA, pval))
  dimnames(table) <- list(
    models,
    c("npar", "Deviance", "Df", "Chisq", "Pr(>Chisq)")
  )
  rownames(table) <- family
  structure(
    table,
    heading = c("Analysis of Deviance Table\n"),
    class = c("anova", "data.frame")
  )
  return(table)
  # }
}


#' Goodness-of-fit diagnostics
#'
#' Warning: EXPERIMENTAL
#' Compute the Kolmogorov-Smirnov
#' test statistic and compare it with a simulated null
#' distribution obtained via a parametric bootstrap.
#'
#' @note The bootstrap scheme requires simulating new data,
#' fitting a parametric model and estimating the nonparametric
#' maximum likelihood estimate for each new sample.
#' This is computationally intensive in large samples.
#'
#' @inheritParams fit_elife
#' @param B number of bootstrap simulations
#' @return a list with elements
#' \itemize{
#' \item \code{stat} the value of the test statistic
#' \item \code{pval} p-value obtained via simulation
#' }
#' @keywords internal
ks_test <- function(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  ltrunc = NULL,
  rtrunc = NULL,
  type = c("right", "left", "interval", "interval2"),
  family = c(
    "exp",
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
    "beardmake"
  ),
  B = 999L,
  arguments = NULL,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = "ks_test",
      call = call,
      arguments = arguments
    )
    return(do.call(ks_test, args = arguments))
  }
  # Exclude doubly interval truncated data
  if (
    isTRUE(all(
      is.matrix(ltrunc),
      is.matrix(rtrunc),
      ncol(ltrunc) == ncol(rtrunc),
      ncol(rtrunc) == 2L
    ))
  ) {
    stop("Doubly interval truncated data not supported")
  }

  family <- match.arg(family)
  type <- match.arg(type)
  n <- length(time)
  survout <- .check_surv(time = time, time2 = time2, event = event, type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  if (thresh[1] > 0) {
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2, FALSE, time > thresh[1])
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    time2 <- time2[ind] - thresh[1]
    status <- status[ind]
    if (!is.null(ltrunc)) {
      #both ltrc and ltrt
      stopifnot(
        "`ltrunc` must be of the same length as `time`." = length(ltrunc) == n
      )
      ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
    }
    if (!is.null(rtrunc)) {
      stopifnot(
        "`rtrunc` must be of the same length as `time`." = length(rtrunc) == n
      )
      rtrunc <- rtrunc[ind] - thresh[1]
    }
  }
  thresh <- 0
  # Fit parametric model
  F0 <- try(fit_elife(
    time = time,
    time2 = time2,
    status = status,
    thresh = thresh,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    type = type,
    family = family
  ))
  if (is.character(F0) || !F0$convergence) {
    stop("Could not estimate the parametric model.")
  }
  # Fit NPMLE of ECDF
  Fn <- npsurv(
    time = time,
    type = "interval",
    event = status,
    ltrunc = ltrunc,
    rtrunc = rtrunc
  )
  # Compute test statistic
  if (family == "gompmake") {
    scale <- c(F0$par[1], F0$par[3])
    shape <- F0$par[2]
  } else {
    scale <- F0$par[1]
    shape <- F0$par[-1]
  }
  ks <- max(abs(
    Fn$cdf(Fn$x) -
      pelife(q = Fn$x, scale = scale, shape = shape, family = family)
  ))
  stat <- rep(NA, B + 1L)
  stat[B + 1] <- ks
  if (is.null(ltrunc) && is.null(rtrunc) && isTRUE(all(status == 1L))) {
    type2 <- "none"
  } else if (
    isTRUE(all(status == 1L)) && (!is.null(rtrunc) || !is.null(ltrunc))
  ) {
    type2 <- "ltrt"
  } else if (is.null(rtrunc) && isTRUE(all(status == c(0L, 1L)))) {
    type2 <- "ltrc"
  }
  for (b in 1:B) {
    bootconv <- FALSE
    while (!bootconv) {
      if (type2 == "ltrt") {
        boottime <- samp_elife(
          n = length(time),
          scale = scale,
          shape = shape,
          lower = ltrunc,
          upper = rtrunc,
          family = family,
          type2 = type2
        )
        bootstatus <- rep(1L, length(time))
      } else if (type2 == "none") {
        boottime <- relife(
          n = length(time),
          scale = scale,
          shape = shape,
          family = family
        )
        bootstatus <- rep(1L, length(time))
      } else if (type2 == "ltrc") {
        boot_sim <- samp_elife(
          n = length(time),
          scale = scale,
          shape = shape,
          lower = ltrunc,
          upper = rtrunc,
          family = family,
          type2 = type2
        )
        boottime <- boot_sim$dat
        boottime2 <- ifelse(boot_sim$rcens, NA, boot_sim$dat)
        bootstatus <- as.integer(!boot_sim$rcens)
      }
      # bootstrap loop
      F0_b <- try(fit_elife(
        time = boottime,
        time2 = boottime2,
        status = bootstatus,
        thresh = thresh,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        type = "right",
        family = family
      ))
      # Fit NPMLE of ECDF
      if (type2 == "ltrc") {
        Fn_b <- try(npsurv(
          time = boottime,
          event = bootstatus,
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          type = "interval"
        ))
      } else {
        Fn_b <- try(npsurv(time = boottime, ltrunc = ltrunc, rtrunc = rtrunc))
      }
      # TODO replace this with function that takes par and returns a list with arguments
      # rate, scale and shape
      # TODO check that all calls to *elife, etc. also have a rate parameter
      if (family == "gompmake") {
        scale_boot <- c(F0_b$par[1], F0_b$par[3])
        shape_boot <- F0_b$par[2]
      } else {
        scale_boot <- F0_b$par[1]
        shape_boot <- F0_b$par[-1]
      }
      # Compute test statistic
      ks_boot <- try(max(abs(
        Fn_b$cdf(Fn_b$x) -
          pelife(
            q = Fn_b$x,
            scale = scale_boot,
            shape = shape_boot,
            family = family
          )
      )))
      if (is.numeric(ks_boot)) {
        stat[b] <- ks_boot
        bootconv <- TRUE
      }
    }
  }
  list(
    stat = ks,
    # null = stat,
    pval = mean(stat >= ks, na.rm = TRUE)
  )
}
