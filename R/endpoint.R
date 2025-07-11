#' Profile likelihood for the endpoint of the generalized Pareto distribution
#'
#' This function returns the profile log likelihood over a grid of values of \code{psi}, the endpoints.
#'
#' @export
#' @inheritParams hazard_elife
#' @param psi mandatory vector of endpoints at which to compute the profile
#' @param confint logical; if \code{TRUE}, return a \code{level} confidence interval instead of a list with the profile log-likelihood components
#' @param level numeric; the level for the confidence intervals
#' @param arguments a named list specifying default arguments of the function that are common to all \code{elife} calls
#' @param ... additional parameters, currently ignored
#' @return a list with the maximum likelihood estimate of the endpoint and the profile log-likelihood
#' @examples
#' set.seed(2023)
#' time <- relife(n = 100, scale = 3, shape = -0.3, family = "gp")
#' endpt <- endpoint.profile(
#'   time = time,
#'   psi = seq(max(time) + 1e-4, max(time) + 40, length.out = 51L))
#' print(endpt)
#' plot(endpt)
#' confint(endpt)
endpoint.profile <- function(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  type = c("right", "left", "interval", "interval2"),
  ltrunc = NULL,
  rtrunc = NULL,
  weights = rep(1, length(time)),
  psi = NULL,
  confint = FALSE,
  level = 0.95,
  arguments = NULL,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = endpoint.profile,
      call = call,
      arguments = arguments
    )
    return(do.call(endpoint.profile, args = arguments))
  }
  stopifnot(
    "Argument \"psi\" must be provided (currently NULL)" = !is.null(psi),
    "Threshold must be positive" = isTRUE(all(thresh >= 0)) &
      length(thresh) == 1,
    "Endpoints must be positive and larger than the threshold" = isTRUE(all(
      psi > thresh
    ))
  )
  psi <- psi - thresh[1]
  type <- match.arg(type)
  # Compute the maximum log-likelihood
  mle <- fit_elife(
    time = time,
    time2 = time2,
    event = event,
    thresh = thresh,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    type = type,
    family = "gp",
    weights = weights
  )
  np <- length(psi)
  param <- matrix(nrow = np, ncol = 2L)
  param[, 1] <- psi
  colnames(param) <- c("endpt", "shape")
  ll <- vector(mode = "numeric", np)
  opt_fun <- function(xi, endpoint) {
    sigma <- -endpoint * xi
    nll_elife(
      par = c(sigma, xi),
      time = time,
      time2 = time2,
      event = event,
      weights = weights,
      thresh = thresh,
      type = type,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp"
    )
  }
  for (i in seq_along(psi)) {
    opt <- optim(
      fn = opt_fun,
      method = "Brent",
      par = -0.1,
      lower = -1,
      upper = -1e-8,
      control = list(reltol = 1e-12),
      endpoint = psi[i]
    )
    ll[i] <- -opt$value
    param[i, 2] <- opt$par
  }
  mle_endpt <- as.numeric(
    thresh + ifelse(mle$par[2] >= 0, Inf, -mle$par[1] / mle$par[2])
  )
  prof <- structure(
    list(
      psi = psi + thresh,
      lambda = param[, 2],
      pll = ll,
      maxpll = mle$loglik,
      mle = mle_endpt,
      psi.max = mle_endpt,
      nexc = mle$nexc,
      param = "endpoint"
    ),
    class = "elife_profile"
  )
  if (is.infinite(mle_endpt)) {
    warning("Maximum likelihood estimate is infinite.")
  }
  return(prof)
}

#' @export
confint.elife_profile <-
  function(object, parm, level = 0.95, ...) {
    confint <- suppressWarnings(try(conf_interv(object, level = level)))
    if (inherits(confint, "try-error")) {
      stop("Could not compute confidence interval")
    } else {
      return(confint)
    }
  }

#' @export
print.elife_profile <- function(x, ...) {
  cat("Parameter:", x$param, "\n")
  cat("Maximum likelihood estimator: ", round(x$psi.max, 3), "\n")
}

#' Threshold stability plots for endpoint
#'
#' This function calls \code{endpoint.profile} to compute
#' the endpoint with profile likelihood confidence intervals
#' at each threshold. It then returns a data frame and
#' a plot
#' @inheritParams endpoint.profile
#' @param plot logical; if \code{TRUE}, return a plot
#' @return a data frame with the threshold, endpoint estimates and profile confidence intervals
#' @export
endpoint.tstab <- function(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  type = c("right", "left", "interval", "interval2"),
  ltrunc = NULL,
  rtrunc = NULL,
  weights = rep(1, length(time)),
  psi = NULL,
  confint = FALSE,
  level = 0.95,
  arguments = NULL,
  plot = TRUE,
  ...
) {
  if (!is.null(arguments)) {
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(
      func = endpoint.tstab,
      call = call,
      arguments = arguments
    )
    return(do.call(endpoint.tstab, args = arguments))
  }
  if (length(thresh) <= 1) {
    stop("Vector \"thresh\" must contain multiple thresholds.")
  }
  profiles <- lapply(thresh, function(thresholds) {
    do.call(
      endpoint.profile,
      args = list(
        time = time,
        time2 = time2,
        event = event,
        type = type,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        weights = weights,
        psi = psi,
        confint = TRUE,
        level = level,
        arguments = arguments,
        thresh = thresholds
      )
    )
  })
  # Compute corresponding confidence intervals
  conf <- data.frame(
    threshold = thresh,
    t(sapply(profiles, confint, level = level))
  )
  class(conf) <- c("elife_tstab_endpoint", "data.frame")
  if (isTRUE(plot)) {
    plot(conf, ...)
  }
  return(invisible(conf))
}

#' Plot profile of endpoint
#'
#' @param x an object of class \code{elife_profile} containing information about the profile likelihood, maximum likelihood and grid of values for the endpoint
#' @param plot.type string indicating whether to use base R for plots or \code{ggplot2}
#' @param plot logical; if \code{TRUE}, creates a plot when \code{plot.type="ggplot"}. Useful for returning \code{ggplot} objects without printing the graphs
#' @param ... additional arguments to pass to \code{plot}, currently ignored
#' @return base R or ggplot object for a plot of the profile log likelihood of the endpoint of the generalized Pareto distribution
#' @export
plot.elife_profile <-
  function(x, plot.type = c("base", "ggplot"), plot = TRUE, ...) {
    plot.type <- match.arg(plot.type)
    plot.type <- match.arg(plot.type)
    if (plot.type == "ggplot") {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      } else {
        warning(
          "`ggplot2` package is not installed. Switching to base R plots."
        )
        plot.type <- "base"
      }
    }
    ind <- which(is.finite(x$pll))
    stopifnot("Maximum log likelihood is not finite" = is.finite(x$maxpll))
    if (plot.type == "base") {
      plot(
        x = x$psi[ind],
        y = x$pll[ind] - x$maxpll,
        type = "l",
        yaxs = "i",
        bty = "l",
        ylim = c(-4, 0.01),
        xlab = x$param,
        ylab = "profile log likelihood"
      )
      abline(h = -qchisq(c(0.95, 0.99), df = 1) / 2, lty = 2)
    } else if (plot.type == "ggplot") {
      g <- ggplot2::ggplot(
        data = data.frame(y = x$pll[ind] - x$maxpll, x = x$psi[ind]),
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]])
      ) +
        ggplot2::geom_hline(
          yintercept = -qchisq(c(0.95, 0.99), df = 1) / 2,
          col = "gray"
        ) +
        ggplot2::geom_line() +
        ggplot2::labs(x = x$param, y = "profile log likelihood") +
        ggplot2::theme_classic()
      if (isTRUE(plot)) {
        print(g)
      }
      return(invisible(g))
    }
  }
#' Plot threshold stability plot for endpoint
#'
#' @param x an object of class \code{elife_tstab_endpoint} containing information about the profile likelihood, maximum likelihood and grid of values for the endpoint
#' @param plot.type string indicating whether to use base R for plots or \code{ggplot2}
#' @param plot logical; if \code{TRUE}, creates a plot when \code{plot.type="ggplot"}. Useful for returning \code{ggplot} objects without printing the graphs
#' @param ... additional arguments to pass to \code{plot}, currently ignored
#' @return base R or ggplot object for a plot of the profile log likelihood of the endpoint of the generalized Pareto distribution
#' @export
plot.elife_tstab_endpoint <-
  function(
    x,
    plot.type = c("base", "ggplot"),
    plot = TRUE,
    ...
  ) {
    args <- list(...)
    plot.type <- match.arg(plot.type)
    plot.type <- match.arg(plot.type)
    if (plot.type == "ggplot") {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      } else {
        warning(
          "`ggplot2` package is not installed. Switching to base R plots."
        )
        plot.type <- "base"
      }
    }
    if (plot.type == "base") {
      plot(
        x = x$threshold,
        y = x$Estimate,
        type = "l",
        yaxs = "i",
        bty = "l",
        ylim = c(min(x$Lower.CI), max(x$Upper.CI)),
        xlab = ifelse(is.null(args$ylab), "threshold", args$ylab),
        ylab = ifelse(is.null(args$ylab), "endpoint", args$ylab)
      )
      lines(x = x$threshold, y = x$Lower.CI, lty = 2)
      lines(x = x$threshold, y = x$Upper.CI, lty = 2)
    } else if (plot.type == "ggplot") {
      g <- ggplot2::ggplot(
        data = as.data.frame(x),
        mapping = ggplot2::aes(
          x = .data[["threshold"]],
          y = .data[["Estimate"]],
          ymin = .data[["Lower.CI"]],
          ymax = .data[["Upper.CI"]]
        )
      ) +
        ggplot2::geom_pointrange() +
        ggplot2::labs(
          x = ifelse(is.null(args$ylab), "threshold", args$ylab),
          y = ifelse(is.null(args$ylab), "endpoint", args$ylab)
        ) +
        ggplot2::theme_classic()
      if (isTRUE(plot)) {
        print(g)
      }
      return(invisible(g))
    }
  }

#' Profile likelihood for generalized Pareto endpoint
#' @export
#' @rdname endpoint.tstab
#' @keywords internal
prof_gp_endpt <- function(
  time,
  time2 = NULL,
  event = NULL,
  thresh = 0,
  type = c("right", "left", "interval", "interval2"),
  ltrunc = NULL,
  rtrunc = NULL,
  weights = rep(1, length(time)),
  psi = NULL,
  confint = FALSE,
  level = 0.95,
  arguments = NULL,
  ...
) {
  .Deprecated(new = "endpoint.profile")
  endpoint.profile(
    time = time,
    time2 = time2,
    event = event,
    thresh = thresh,
    type = type,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    weights = weights,
    psi = psi,
    confint = confint,
    level = level,
    arguments = arguments,
    ...
  )
}
