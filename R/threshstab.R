#' Threshold stability plots
#'
#' The generalized Pareto and exponential distribution
#' are threshold stable. This property, which is used
#' for extrapolation purposes, can also be used to diagnose
#' goodness-of-fit: we expect the parameters \eqn{\xi} and \eqn{\tilde{\sigma} = \sigma + \xi u}
#' to be constant over a range of thresholds. The threshold stability
#' plot consists in plotting maximum likelihood estimates with pointwise confidence interval.
#' This function handles interval truncation and right-censoring.
#'
#' @details The shape estimates are constrained
#' @inheritParams nll_elife
#' @param family string; distribution, either generalized Pareto (\code{gp}) or exponential (\code{exp})
#' @param method string; the type of pointwise confidence interval, either Wald (\code{wald}) or profile likelihood (\code{profile})
#' @param level probability level for the pointwise confidence intervals
#' @param plot.type string; either \code{base} for base R plots or \code{ggplot} for \code{ggplot2} plots
#' @param which.plot string; which parameters to plot;
#' @param plot logical; should a plot be returned alongside with the estimates? Default to \code{TRUE}
#' @seealso \code{tstab.gpd} from package \code{mev}, \code{gpd.fitrange} from package \code{ismev} or \code{tcplot} from package \code{evd}, among others.
#' @importFrom utils packageVersion
#' @export
#' @return an invisible list with pointwise estimates and confidence intervals for the scale and shape parameters
#' @examples
#' set.seed(1234)
#' n <- 100L
#' x <- samp_elife(n = n,
#'                 scale = 2,
#'                 shape = -0.2,
#'                 lower = low <- runif(n),
#'                 upper = upp <- runif(n, min = 3, max = 20),
#'                 type2 = "ltrt",
#'                 family = "gp")
#' tstab_plot <- tstab(time = x,
#'                     ltrunc = low,
#'                    rtrunc = upp,
#'                    thresh = quantile(x, seq(0, 0.5, length.out = 4)))
#' plot(tstab_plot, plot.type = "ggplot")
tstab <- function(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  ltrunc = NULL,
  rtrunc = NULL,
  type = c("right", "left", "interval", "interval2"),
  family = c('gp', 'exp'),
  method = c("wald", "profile"),
  level = 0.95,
  plot = TRUE,
  plot.type = c("base", "ggplot"),
  which.plot = c("scale", "shape"),
  weights = NULL,
  arguments = NULL,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = tstab,
      call = call,
      arguments = arguments
    )
    return(do.call(tstab, args = arguments))
  }
  family <- match.arg(family)
  method <- match.arg(method)
  type <- match.arg(type)
  if (plot) {
    plot.type <- match.arg(plot.type)
    which.plot <- match.arg(
      which.plot,
      choices = c("scale", "shape"),
      several.ok = TRUE
    )
  }
  wald_confint <- function(par, std.error, level = 0.95) {
    c(par[1], par[1] + qnorm(0.5 + level / 2) * std.error[1] * c(-1, 1))
  }
  stopifnot("Provide multiple thresholds" = length(thresh) > 1)
  scale_par_mat <- matrix(NA, nrow = length(thresh), ncol = 3)
  if (family == "gp") {
    shape_par_mat <- matrix(NA, nrow = length(thresh), ncol = 3)
  }
  for (i in seq_along(thresh)) {
    opt_mle <- fit_elife(
      time = time,
      time2 = time2,
      event = event,
      thresh = thresh[i],
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      type = type,
      family = family,
      export = TRUE,
      weights = weights
    )
    opt_mle$mdat <- max(c(opt_mle$time, opt_mle$time2), na.rm = TRUE)
    if (method == "profile") {
      if (family == "gp") {
        if ("shape" %in% which.plot) {
          shape_i <- try(prof_gp_shape(
            mle = opt_mle,
            time = time,
            time2 = time2,
            event = event,
            thresh = thresh[i],
            ltrunc = ltrunc,
            rtrunc = rtrunc,
            type = type,
            level = level,
            weights = weights
          ))
          if (!inherits(shape_i, "try.error")) {
            shape_par_mat[i, ] <- shape_i
          }
        }
        if ("scale" %in% which.plot) {
          scalet_i <- try(prof_gp_scalet(
            mle = opt_mle,
            time = time,
            time2 = time2,
            event = event,
            thresh = thresh[i],
            ltrunc = ltrunc,
            rtrunc = rtrunc,
            type = type,
            level = level,
            weights = weights
          ))
          if (!inherits(scalet_i, "try-error")) {
            scale_par_mat[i, ] <- scalet_i
          }
        }
      } else if (family == "exp") {
        scale_i <- try(prof_exp_scale(
          mle = opt_mle,
          time = time,
          time2 = time2,
          event = event,
          thresh = thresh[i],
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          type = type,
          level = level
        ))
        if (!inherits(scale_i, "try-error")) {
          scale_par_mat[i, ] <- scale_i
        }
      }
    } else if (method == "wald") {
      if (family == "gp") {
        sigma_t_mle <- opt_mle$par[1] - opt_mle$par[2] * thresh[i]
        dV <- matrix(c(1, -thresh[i]), ncol = 1)
        if (!is.null(opt_mle$vcov)) {
          v <- t(dV) %*% opt_mle$vcov %*% dV
          shape_par_mat[i, ] <- wald_confint(
            par = opt_mle$par[2],
            opt_mle$std.error[2],
            level = level
          )
          scale_par_mat[i, ] <- wald_confint(
            par = sigma_t_mle,
            std.error = sqrt(v),
            level = level
          )
        }
      } else if (family == "exp") {
        if (!is.null(opt_mle$vcov)) {
          scale_par_mat[i, ] <- wald_confint(
            par = opt_mle$par,
            opt_mle$std.error
          )
        }
      }
    }
  }
  colnames(scale_par_mat) <- c("estimate", "lower", "upper")
  if (family == "gp") {
    colnames(shape_par_mat) <- c("estimate", "lower", "upper")
    if (!"scale" %in% which.plot) {
      scale_par_mat <- NULL
    }
    if (!"shape" %in% which.plot) {
      shape_par_mat <- NULL
    }
    res <- structure(
      list(
        family = family,
        scale = scale_par_mat,
        shape = shape_par_mat,
        thresh = thresh
      ),
      class = "elife_tstab"
    )
  } else if (family == "exp") {
    res <- structure(
      list(family = family, scale = scale_par_mat, thresh = thresh),
      class = "elife_tstab"
    )
  }
  if (plot) {
    plot(res, plot.type = plot.type, which.plot = which.plot, plot = TRUE)
  }
  #TODO check this for left-truncated data
  res$nexc <- vapply(
    thresh,
    function(u) {
      sum(time > u, na.rm = TRUE)
    },
    integer(1)
  )
  invisible(res)
}

#' Threshold stability plots
#' @param x an object of class \code{elife_tstab} containing
#' the fitted parameters as a function of threshold
#' @inheritParams plot.elife_par
#' @export
plot.elife_tstab <- function(
  x,
  plot.type = c("base", "ggplot"),
  which.plot = c("scale", "shape"),
  plot = TRUE,
  ...
) {
  object <- x
  plot.type <- match.arg(plot.type)
  if (plot.type == "ggplot") {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
    } else {
      warning("`ggplot2` package is not installed. Switching to base R plots.")
      plot.type <- "base"
    }
  }
  which.plot <- match.arg(
    which.plot,
    choices = c("scale", "shape"),
    several.ok = TRUE
  )
  if (plot.type == "ggplot") {
    ggplot_thstab <- function(x, thresh, ylab) {
      stopifnot(all.equal(colnames(x), c("estimate", "lower", "upper")))
      g <- ggplot2::ggplot(
        data = as.data.frame(cbind(thresh = thresh, x)),
        ggplot2::aes(x = .data[["thresh"]], y = .data[["estimate"]])
      ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
          size = 0.5,
          shape = 20
        ) +
        ggplot2::labs(
          x = "threshold",
          y = ylab,
          main = "threshold stability plot"
        ) +
        ggplot2::theme_classic() #+
      # scale_x_continuous(breaks = thresh, minor_breaks = NULL)
      return(invisible(g))
    }
    graphs <- list()
    if (object$family == "gp") {
      if ("scale" %in% which.plot) {
        g1 <- ggplot_thstab(
          x = object$scale,
          thresh = object$thresh,
          ylab = "modified scale"
        )
        if (length(which.plot) == 1L && plot) {
          get(
            ifelse(
              packageVersion("ggplot2") >= "3.5.2.9001",
              "print.ggplot2::ggplot",
              "print.ggplot"
            ),
            envir = loadNamespace("ggplot2")
          )(g1)
        }
        graphs$g1 <- g1
      }
      if ("shape" %in% which.plot) {
        g2 <- ggplot_thstab(
          x = object$shape,
          thresh = object$thresh,
          ylab = "shape"
        )
        if (length(which.plot) == 1L && plot) {
          get(
            ifelse(
              packageVersion("ggplot2") >= "3.5.2.9001",
              "print.ggplot2::ggplot",
              "print.ggplot"
            ),
            envir = loadNamespace("ggplot2")
          )(g2)
        }
        graphs$g2 <- g2
      }
      if (length(which.plot) == 2L && plot) {
        lapply(
          list(g1, g2),
          get(
            ifelse(
              packageVersion("ggplot2") >= "3.5.2.9001",
              "print.ggplot2::ggplot",
              "print.ggplot"
            ),
            envir = loadNamespace("ggplot2")
          )
        )
      }
    } else if (object$family == "exp") {
      g1 <- ggplot_thstab(
        x = object$scale,
        thresh = object$thresh,
        ylab = "scale"
      )
      if (plot) {
        get(
          ifelse(
            packageVersion("ggplot2") >= "3.5.2.9001",
            "print.ggplot2::ggplot",
            "print.ggplot"
          ),
          envir = loadNamespace("ggplot2")
        )(g1)
      }
      graphs$g1 <- g1
    }
    return(invisible(graphs))
  } else {
    # Default to base plot
    base_tstab_plot <- function(x, thresh, ylab = "scale") {
      rangex <- range(x)
      args <- list(
        type = "p",
        pch = 19,
        bty = "l",
        ylab = ylab,
        xlab = "threshold",
        ylim = rangex
      )
      ellipsis <- list(...)
      args[names(ellipsis)] <- ellipsis
      args$x <- thresh
      args$y <- x[, 1]
      do.call("plot.default", args = args)
      for (i in seq_along(thresh)) {
        arrows(
          x0 = thresh[i],
          y0 = x[i, 2],
          y1 = x[i, 3],
          length = 0.02,
          angle = 90,
          code = 3,
          lwd = 2
        )
      }
    }
    if (object$family == "exp") {
      base_tstab_plot(x = object$scale, thresh = object$thresh)
    } else {
      #generalized Pareto
      if ("scale" %in% which.plot) {
        base_tstab_plot(
          x = object$scale,
          thresh = object$thresh,
          ylab = "modified scale"
        )
      }
      if ("shape" %in% which.plot) {
        base_tstab_plot(
          x = object$shape,
          thresh = object$thresh,
          ylab = "shape"
        )
      }
    }
  }
}


#' Profile log likelihood for the shape parameter of the generalized Pareto distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param confint logical; if \code{TRUE} (default), return confidence intervals rather than list
#' @param level level of the confidence interval
#' @keywords internal
#' @export
#' @return if \code{confint=TRUE}, a vector of length three with the maximum likelihood of the shape and profile-based confidence interval
prof_gp_shape <-
  function(
    mle = NULL,
    time,
    time2 = NULL,
    event = NULL,
    thresh,
    ltrunc = NULL,
    rtrunc = NULL,
    type = c("right", "left", "interval", "interval2"),
    level = 0.95,
    psi = NULL,
    weights = NULL,
    confint = TRUE,
    arguments = NULL,
    ...
  ) {
    if (!is.null(arguments)) {
      call <- match.call(expand.dots = FALSE)
      arguments <- check_arguments(
        func = prof_gp_shape,
        call = call,
        arguments = arguments
      )
      return(do.call(prof_gp_shape, args = arguments))
    }
    if (is.null(weights)) {
      weights <- rep(1, length(time))
    }
    type <- match.arg(type)
    stopifnot(
      "Provide a single value for the level" = length(level) == 1L,
      "Level should be a probability" = level < 1 && level > 0,
      "Provide a single threshold" = length(thresh) == 1L
    )
    if (!is.null(mle)) {
      stopifnot(
        "`mle` should be an object of class `elife_par` as returned by `fit_elife`" = inherits(
          mle,
          "elife_par"
        )
      )
    } else {
      mle <- fit_elife(
        time = time,
        time2 = time2,
        event = event,
        thresh = thresh,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        type = type,
        family = "gp",
        export = TRUE,
        weights = weights
      )
      mle$mdat <- max(c(mle$time, mle$time2), na.rm = TRUE)
    }
    if (is.null(psi)) {
      psi <- seq(-0.99, 2, length = 100L)
    } else {
      stopifnot(
        "`psi` should be a numeric vector" = is.numeric(psi),
        "Grid of values for `psi` do not include the maximum likelihood estimate." = min(
          psi
        ) <
          mle$par[2] &
          max(psi) > mle$par[2]
      )
    }
    mdat <- mle$mdat
    dev <- vapply(
      psi,
      function(xi) {
        opt <- optimize(
          f = function(lambda) {
            nll_elife(
              par = c(lambda, xi),
              time = time,
              time2 = time2,
              event = event,
              thresh = thresh,
              type = type,
              ltrunc = ltrunc,
              rtrunc = rtrunc,
              family = "gp",
              weights = weights
            )
          },
          interval = c(ifelse(xi < 0, mdat * abs(xi), 1e-8), 10 * mdat),
          tol = 1e-10
        )
        c(-2 * opt$objective, opt$minimum)
      },
      FUN.VALUE = numeric(2)
    )
    prof <- list(
      psi = psi,
      lambda = dev[2, ],
      pll = -2 * mle$loglik + dev[1, ],
      maxpll = 0,
      mle = mle$par,
      psi.max = mle$par['shape'],
      std.error = sqrt(mle$vcov[2, 2])
    )
    if (confint) {
      conf_interv(prof, level = level, print = FALSE)
    } else {
      prof
    }
  }

#' Profile log likelihood for the transformed scale parameter of the generalized Pareto distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param confint logical; if \code{TRUE} (default), return confidence intervals rather than list
#' @param level level of the confidence interval
#' @keywords internal
#' @export
#' @return a confidence interval or a list with profile values
prof_gp_scalet <-
  function(
    mle = NULL,
    time,
    time2 = NULL,
    event = NULL,
    thresh = 0,
    ltrunc = NULL,
    rtrunc = NULL,
    type = c("right", "left", "interval", "interval2"),
    level = 0.95,
    psi = NULL,
    weights = NULL,
    confint = TRUE,
    arguments = NULL,
    ...
  ) {
    if (!is.null(arguments)) {
      call <- match.call(expand.dots = FALSE)
      arguments <- check_arguments(
        func = prof_gp_scalet,
        call = call,
        arguments = arguments
      )
      return(do.call(prof_gp_scalet, args = arguments))
    }

    if (is.null(weights)) {
      weights <- rep(1, length(time))
    }
    type <- match.arg(type)
    stopifnot(
      "Provide a single value for the level" = length(level) == 1L,
      "Level should be a probability" = level < 1 && level > 0,
      "Provide a single threshold" = length(thresh) == 1L
    )
    if (!is.null(mle)) {
      stopifnot(
        "`mle` should be an object of class `elife_par` as returned by `fit_elife`" = inherits(
          mle,
          "elife_par"
        )
      )
    } else {
      mle <- fit_elife(
        time = time,
        time2 = time2,
        event = event,
        thresh = thresh,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        type = type,
        family = "gp",
        export = TRUE,
        weights = weights
      )
      mle$mdat <- max(c(mle$time, mle$time2), na.rm = TRUE)
    }
    mdat <- mle$mdat
    sigma_t_mle <- mle$par[1] - mle$par[2] * thresh
    dV <- matrix(c(1, -thresh), ncol = 1)
    sigma_t_se <- sqrt(as.numeric(t(dV) %*% mle$vcov %*% dV))
    if (is.null(psi)) {
      psi <- sigma_t_mle +
        seq(-5 * sigma_t_se, 10 * sigma_t_se, length.out = 101)
    } else {
      stopifnot(
        "`psi` should be a numeric vector" = is.numeric(psi),
        "Grid of values for `psi` do not include the maximum likelihood estimate." = min(
          psi
        ) <
          sigma_t_mle &
          max(psi) > sigma_t_mle
      )
    }
    psi <- psi[psi > 0]
    # Optimize is twice as fast as Rsolnp...
    dev <- t(vapply(
      psi,
      function(scalet) {
        opt <- optimize(
          f = function(xi) {
            nll_elife(
              par = c(scalet + xi * thresh, xi),
              time = time,
              time2 = time2,
              event = event,
              thresh = thresh,
              type = type,
              ltrunc = ltrunc,
              rtrunc = rtrunc,
              family = "gp",
              weights = weights
            )
          },
          interval = c(
            pmax(-1, ifelse(thresh == 0, -scalet / thresh, -scalet / mdat)),
            3
          ),
          tol = 1e-10
        )
        c(-2 * opt$objective, opt$minimum)
      },
      FUN.VALUE = numeric(2)
    ))
    prof <- list(
      psi = psi,
      pll = -2 * mle$loglik + dev[, 1],
      maxpll = 0,
      mle = mle$par,
      psi.max = sigma_t_mle,
      std.error = sigma_t_se
    )
    if (confint) {
      conf_interv(prof, level = level, print = FALSE)
    } else {
      prof
    }
  }


#' Profile log likelihood for the scale parameter of the exponential distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param confint logical; if \code{TRUE} (default), return confidence intervals rather than list
#' @param level level of the confidence interval
#' @export
#' @return a vector of length three with the maximum likelihood of the scale and profile-based confidence interval
#' @keywords internal
prof_exp_scale <- function(
  mle = NULL,
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  ltrunc = NULL,
  rtrunc = NULL,
  type = c("right", "left", "interval", "interval2"),
  level = 0.95,
  psi = NULL,
  weights = NULL,
  confint = TRUE,
  arguments = NULL,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = prof_exp_scale,
      call = call,
      arguments = arguments
    )
    return(do.call(prof_exp_scale, args = arguments))
  }
  type <- match.arg(type)
  stopifnot(
    "Provide a single value for the level" = length(level) == 1L,
    "Level should be a probability" = level < 1 && level > 0,
    "Provide a single threshold" = length(thresh) == 1L
  )
  if (is.null(weights)) {
    weights <- rep(1, length(time))
  }
  if (!is.null(mle)) {
    stopifnot(
      "`mle` should be an object of class `elife_par` as returned by `fit_elife`" = inherits(
        mle,
        "elife_par"
      )
    )
  } else {
    mle <- fit_elife(
      time = time,
      time2 = time2,
      event = event,
      thresh = thresh,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      type = type,
      family = "exp",
      weights = weights
    )
  }
  if (is.null(psi)) {
    psi <- mle$par +
      seq(
        pmax(-mle$par + 1e-4, -4 * mle$std.error),
        4 * mle$std.error,
        length.out = 201
      )
  } else {
    stopifnot(
      "`psi` should be a numeric vector" = is.numeric(psi),
      "Grid of values for `psi` do not include the maximum likelihood estimate." = min(
        psi
      ) <
        mle$par &
        max(psi) > mle$par
    )
  }

  pll <- vapply(
    psi,
    function(scale) {
      nll_elife(
        par = scale,
        time = time,
        time2 = time2,
        event = event,
        thresh = thresh,
        type = type,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        family = "exp",
        weights = weights
      )
    },
    FUN.VALUE = numeric(1)
  )
  prof <- list(
    psi = psi,
    pll = -2 * (mle$loglik + pll),
    maxpll = 0,
    psi.max = mle$par,
    std.error = mle$std.error,
    mle = mle$par
  )
  if (confint) {
    conf_interv(prof, level = level)
  } else {
    prof
  }
}


#' Confidence intervals for profile likelihoods
#'
#' This code is adapted from the \code{mev} package.
#' @param object a list containing information about the profile likelihood in the same format as the \code{hoa} package
#' @param level probability level of the confidence interval
#' @param prob vector of length 2 containing the bounds, by default double-sided
#' @param print logical indicating whether the intervals are printed to the console
#' @param ... additional arguments passed to the function
#' @return a table with confidence intervals.
#' @keywords internal
conf_interv <- function(
  object,
  level = 0.95,
  prob = c((1 - level) / 2, 1 - (1 - level) / 2),
  print = FALSE,
  ...
) {
  if (!isTRUE(all.equal(diff(prob), level, check.attributes = FALSE))) {
    warning("Incompatible arguments: `level` does not match `prob`.")
  }
  args <- list(...)
  if ("warn" %in% names(args) && is.logical(args$warn)) {
    warn <- args$warn
  } else {
    warn <- TRUE
  }
  if (length(prob) != 2) {
    stop("`prob` must be a vector of size 2")
    prob <- sort(prob)
  }
  qulev <- qnorm(1 - prob)
  conf <- rep(0, 3)
  if (is.null(object$pll) && is.null(object$r)) {
    stop(
      "Object should contain arguments `pll` or `r` in order to compute confidence intervals."
    )
  }
  if (is.null(object$r)) {
    object$r <- sign(object$psi.max - object$psi) *
      sqrt(2 * (object$maxpll - object$pll))
  } else {
    object$r[is.infinite(object$r)] <- NA
  }
  if (is.null(object$normal)) {
    object$normal <- c(object$psi.max, object$std.error)
  }
  fit.r <- stats::smooth.spline(
    x = na.omit(cbind(object$r, object$psi)),
    cv = FALSE
  )
  pr <- predict(fit.r, c(0, qulev))$y
  pr[1] <- object$normal[1]
  conf <- pr
  if (warn) {
    if (!any(object$r > qnorm(prob[1]))) {
      warning(
        "Extrapolating the lower confidence interval for the profile likelihood ratio test"
      )
    }
    if (!any(object$r < qnorm(prob[2]))) {
      warning(
        "Extrapolating the upper confidence interval for the profile likelihood ratio test"
      )
    }
  }
  names(conf) <- c("Estimate", "Lower CI", "Upper CI")
  conf[2] <- ifelse(conf[2] > conf[1], NA, conf[2])
  conf[3] <- ifelse(conf[3] < conf[1], NA, conf[3])
  if (print) {
    cat("Point estimate for the parameter of interest psi:\n")
    cat("Maximum likelihood          :", round(object$psi.max, 3), "\n")
    cat("\n")
    cat("Confidence intervals, levels :", prob, "\n")
    cat(
      "Wald intervals               :",
      round(
        object$psi.max +
          sort(qulev) * object$std.error,
        3
      ),
      "\n"
    )
    cat("Profile likelihood           :", round(conf[2:3], 3), "\n")
  }
  return(invisible(conf))
}
