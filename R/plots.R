#' Goodness-of-fit plots for parametric models
#'
#' Because of censoring and truncation, the plotting
#' positions must be adjusted accordingly.
#' For right-censored data, the methodology is described
#' in Waller & Turnbull (1992). Only non-censored observations are
#' displayed, which can create distortion.
#'
#' For truncated data, we first estimate the distribution function
#' nonparametrically, \eqn{F_n}. The uniform plotting positions of the data
#' \deqn{v_i = [F_n(y_i) - F_n(a_i)]/[F_n(b_i) - F_n(a_i)].}
#' For probability-probability plots, the empirical quantiles are transformed
#' using the same transformation, with \eqn{F_n} replaced by the postulated or estimated
#' distribution function \eqn{F_0}.
#' For quantile-quantile plots, the plotting positions \eqn{v_i} are mapped back
#' to the data scale viz. \deqn{F_0^{-1}\{F_0(a_i) + v_i[F_0(b_i) - F_0(a_i)]\}}
#' When data are truncated and observations are mapped back to the untruncated scale (with, e.g., \code{exp}), the plotting positions need not be in the same order as the order statistics of the data.
#'
#' @export
#' @param x a parametric model of class \code{elife_par}
#' @param plot.type string, one of \code{base} for base R or \code{ggplot}
#' @param which.plot vector of string indicating the plots, among \code{pp} for probability-probability plot, \code{qq} for quantile-quantile plot, \code{erp} for empirically rescaled plot (only for censored data), \code{exp} for standard exponential quantile-quantile plot or \code{tmd} for Tukey's mean difference plot, which is a variant of the Q-Q plot in which we map the pair \eqn{(x,y)} is mapped to \code{((x+y)/2,y-x)} are detrended
#' @param confint logical; if \code{TRUE}, creates uncertainty diagnostic via a parametric bootstrap
#' @param plot logical; if \code{TRUE}, creates a plot when \code{plot.type="ggplot"}. Useful for returning \code{ggplot} objects without printing the graphs
#' @param ... additional arguments, currently ignored by the function.
#' @importFrom rlang .data
#' @return The function produces graphical goodness-of-fit plots using base R or ggplot objects (returned as an invisible list).
#' @examples
#' set.seed(1234)
#' samp <- samp_elife(
#'  n = 200,
#'  scale = 2,
#'  shape = 0.3,
#'  family = "gomp",
#'  lower = 0, upper = runif(200, 0, 10),
#'  type2 = "ltrc")
#' fitted <- fit_elife(
#'  time = samp$dat,
#'  thresh = 0,
#'  event = ifelse(samp$rcens, 0L, 1L),
#'  type = "right",
#'  family = "exp",
#'  export = TRUE)
#' plot(fitted, plot.type = "ggplot")
#' # Left- and right-truncated data
#' n <- 40L
#' samp <- samp_elife(
#'  n = n,
#'  scale = 2,
#'  shape = 0.3,
#'  family = "gp",
#'  lower = ltrunc <- runif(n),
#'  upper = rtrunc <- ltrunc + runif(n, 0, 15),
#'  type2 = "ltrt")
#' fitted <- fit_elife(
#'  time = samp,
#'  thresh = 0,
#'  ltrunc = ltrunc,
#'  rtrunc = rtrunc,
#'  family = "gp",
#'  export = TRUE)
#' plot(fitted,  which.plot = "tmd")
plot.elife_par <- function(x,
                           plot.type = c("base","ggplot"),
                           which.plot = c("pp","qq"),
                           confint = c("none", "pointwise", "simultaneous"),
                           plot = TRUE, ...){
  plot.type <- match.arg(plot.type)
  confint <- match.arg(confint)
  if(plot.type == "ggplot"){
    if(requireNamespace("ggplot2", quietly = TRUE)){
    } else{
      warning("`ggplot2` package is not installed. Switching to base R plots.")
      plot.type <- "base"
    }
  }
  object <- x
  stopifnot("Object should be of class `elife_par`" = inherits(object, what = "elife_par"))
  if(is.null(object$time)){
    stop("Object created using a call to `fit_elife` should include the data (`export=TRUE`).")
  }
  which.plot <- match.arg(which.plot, choices = c("pp","qq","erp","exp","tmd"), several.ok = TRUE)
  # Fit a nonparametric survival function (Turnbull, 1976)
  if(is.null(object$rtrunc)){
    object$rtrunc <- rep(Inf, length(object$time))
  }
  if(is.null(object$ltrunc)){
    object$ltrunc <- rep(0, length(object$time))
  }
  np <- np_elife(arguments = object,
                 method = "em",
                 thresh = 0)
  if(!np$convergence){
    warning("Nonparametric maximum likelihood routine did not converge.")
  }
  # Create a weighted empirical CDF
  ecdffun <- np$cdf
  dat <- object$time[object$status == 1L]
  cens <- length(dat) != length(object$time)
  if(!is.null(object$ltrunc)){
    if(length(object$ltrunc) == 1L){
      ltrunc <- rep(object$ltrunc, length.out = sum(object$status == 1L))
    } else{
      stopifnot("Lower truncation limit should be a vector of length 1 or n."  = length(object$ltrunc) == length(object$time))
    ltrunc <- object$ltrunc[object$status == 1L]
    }
  } else{
    ltrunc <- NULL
  }
  if(!is.null(object$rtrunc)){
    if(length(object$rtrunc) == 1L){
      rtrunc <- rep(object$rtrunc, length.out = sum(object$status == 1L))
    } else{
      stopifnot("Upper truncation limit should be a vector of length 1 or n."  = length(object$rtrunc) == length(object$time))
      rtrunc <- object$rtrunc[object$status == 1L]
    }
  } else{
    rtrunc <- NULL
  }
  xpos <- length(dat) * ecdffun(dat) / (length(dat) + 1L)
  if(!is.null(ltrunc) && is.null(rtrunc)){
    xpos <- pmax(0,(xpos - ecdffun(ltrunc)))/(1-ecdffun(ltrunc))
  } else if(!is.null(rtrunc) && is.null(ltrunc)){
    xpos <- xpos/ecdffun(rtrunc)
  } else if(!is.null(rtrunc) && !is.null(ltrunc)){
    xpos <- pmax(0,(xpos - ecdffun(ltrunc)))/(ecdffun(rtrunc) - ecdffun(ltrunc))
  }
  scale <- as.numeric(object$par[1])
  if(object$family == "gompmake"){
    shape <- as.numeric(object$par[2])
    scale <- as.numeric(object$par[-2])
  } else{
    scale <- as.numeric(object$par[1])
    shape <- as.numeric(object$par[-1])
  }
  pmod <- function(q, scale, shape, family){
    switch(object$family,
           exp = pexp(q = q, rate = 1/scale),
           gp = pgpd(q = q, loc = 0, scale = scale, shape = shape),
           gomp = pgomp(q = q, scale = scale, shape = shape),
           gompmake = pgompmake(q = q, scale = scale[1], shape = shape, lambda = scale[2]),
           weibull = pweibull(q = q, shape = shape, scale = scale),
           extgp = pextgp(q = q, scale = scale, shape1 = shape[1], shape2 = shape[2]),
           gppiece = pgppiece(q = q, scale = scale, shape = shape, thresh = object$thresh)
    )
  }
  qmod <- function(p, scale, shape, family){
    switch(object$family,
           exp = qexp(p = p, rate = 1/scale),
           gp = qgpd(p = p, loc = 0, scale = scale, shape = shape),
           gomp = qgomp(p = p, scale = scale, shape = shape),
           gompmake = qgompmake(p = p, scale = scale[1], shape = shape, lambda = scale[2]),
           weibull = qweibull(p = p, shape = shape, scale = scale),
           extgp = qextgp(p = p, scale = scale, shape1 = shape[1], shape2 = shape[2]),
           gppiece = qgppiece(p = p, scale = scale, shape = shape, thresh = object$thresh)
    )
  }
  ypos <- pmod(q = dat, scale = scale, shape = shape, family = family)
  if(!is.null(ltrunc) && is.null(rtrunc)){
    ypos <- (ypos - pmod(q = ltrunc, scale = scale, shape = shape, family = family))/(1-pmod(q = ltrunc, scale = scale, shape = shape, family = family))
  } else if(!is.null(rtrunc) && is.null(ltrunc)){
    ypos <- ypos/pmod(q = rtrunc, scale = scale, shape = shape, family = family)
  } else if(!is.null(rtrunc) && !is.null(ltrunc)){
    ypos <- (ypos - pmod(q = ltrunc, scale = scale, shape = shape, family = family))/(pmod(q = rtrunc, scale = scale, shape = shape, family = family) - pmod(q = ltrunc, scale = scale, shape = shape, family = family))
  }
  if(any(c("qq","tmd","erp") %in% which.plot)){
    txpos <- xpos
    if(!is.null(ltrunc) && !is.null(rtrunc)){
      txpos <- xpos*pmod(q = rtrunc, scale = scale, shape = shape, family = object$family) +
        (1-xpos)*pmod(q = ltrunc, scale = scale, shape = shape, family = object$family)
    } else if(is.null(ltrunc) && !is.null(rtrunc)){
      txpos <- xpos*pmod(q = rtrunc, scale = scale, shape = shape, family = object$family)
    } else if(is.null(rtrunc) && !is.null(ltrunc)){
      txpos <- xpos + (1-xpos)*pmod(q = ltrunc, scale = scale, shape = shape, family = object$family)
    }
  }
  if("erp" %in% which.plot){
    if(isTRUE(all(object$status == 1L))){
      warning("`erp` is only useful for censored data. Use `which.plot = pp` for specifying the equivalent plot.")
      if("pp" %in% which.plot){
        which.plot <- which.plot[which.plot == "erp"]
      } else {
        which.plot[which.plot == "erp"] <- "pp"
      }
    } else{
      np2 <- np_elife(time = dat,
                      ltrunc = ltrunc,
                      rtrunc = rtrunc)
      # Create a weighted empirical CDF
      ecdffun2 <- np2$cdf
    }
  }
  if(plot.type == "base"){
    for(pl in which.plot){
      if(pl == "pp"){
        plot(y = ypos,
             x = xpos,
             bty = "l",
             pch = 20,
             xlab = "theoretical quantiles",
             ylab = "empirical quantiles",
             panel.first = {abline(a = 0, b = 1, col = "gray")})
      } else if(pl == "exp"){
        # exponential q-q plot
        plot(y = -log(1-ypos),
             x = -log(1-xpos),
             bty = "l",
             pch = 20,
             xlab = "theoretical quantiles",
             ylab = "empirical quantiles",
             panel.first = {abline(a = 0, b = 1, col = "gray")})
      } else if(pl == "qq"){
        plot(y = dat,
             x = qmod(p = txpos, scale = scale, shape = shape, family = object$family),
             bty = "l",
             pch = 20,
             xlab = "theoretical quantiles",
             ylab = "empirical quantiles",
             panel.first = {abline(a = 0, b = 1, col = "gray")})
      } else if(pl == "tmd"){
        yp <- dat
        xp <- qmod(p = txpos, scale = scale, shape = shape, family = object$family)
        plot(y = yp - xp,
             x = (xp + yp) / 2,
             bty = "l",
             pch = 20,
             xlab = "average quantile",
             ylab = "quantile difference",
             panel.first = {abline(h = 0, col = "gray")})
      } else if(pl == "erp" && cens){
        plot(y = ecdffun2(dat),
             x = ecdffun2(qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
             bty = "l",
             pch = 20,
             xlab = "theoretical quantiles",
             ylab = "empirical quantiles",
             panel.first = {abline(a = 0, b = 1, col = "gray")})
      }
    }
    return(invisible(NULL))
  } else if(plot.type == "ggplot"){
    pl_list <- list()
    for(pl in which.plot){
      if(pl == "pp"){
        pl_list[["pp"]] <-
          ggplot2::ggplot(data = data.frame(y = ypos, x = xpos)) +
          ggplot2::geom_abline(intercept = 0, slope = 1, col = "gray") +
          ggplot2::geom_point(mapping = ggplot2::aes(
            x = .data[["x"]],
            y = .data[["y"]])) +
          ggplot2::labs(x = "theoretical quantiles",
               y = "empirical quantiles") +
          ggplot2::theme_classic()
      } else if(pl == "exp"){
        pl_list[["exp"]] <-
          ggplot2::ggplot(data = data.frame(y = -log(1-ypos),
                                   x = -log(1-xpos)),
                 mapping = ggplot2::aes(x = .data[["x"]],
                                        y = .data[["y"]])) +
          ggplot2::geom_abline(intercept = 0, slope = 1, col = "gray") +
          ggplot2::geom_point() +
          ggplot2::labs(x = "theoretical quantiles",
               y = "empirical quantiles") +
          ggplot2::theme_classic()
      } else if(pl == "qq"){
        # if(confint && object$type != "ltrc"){
        #     # TODO fix this
        #  confint_qq <- uq1_qqplot_elife(B = 1999L,
        #                               n = length(object$time)
        #                                 par = object$par,
        #                                 lower = ltrunc,
        #                                 upper = rtrunc,
        #                                 type2 = object$type,
        #                                 family = object$family)
        # }
        # if(!confint){
        pl_list[["qq"]] <-
          ggplot2::ggplot(data = data.frame(y = dat,
                                   x = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
                 mapping = ggplot2::aes(x = .data[["x"]],
                                        y = .data[["y"]])) +
          ggplot2::geom_abline(intercept = 0, slope = 1, col = "gray") +
          ggplot2::geom_point() +
          ggplot2::labs(x = "theoretical quantiles",
               y = "empirical quantiles") +
          ggplot2::theme_classic()
        # } else if(confint){
        #   pl_list[["qq"]] <-
        #     ggplot(data = data.frame(y = dat,
        #                              x = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
        #            mapping = aes(x = x, y = y)) +
        #     geom_abline(intercept = 0, slope = 1, col = "gray") +
        #     labs(x = "theoretical quantiles",
        #          y = "empirical quantiles") +
        #     geom_line(data = data.frame(x = confint_qq$point[1,][order(object$dat)], y = sort(object$dat)), mapping = aes(x=x,y=y), col = "grey") +
        #     geom_line(data = data.frame(x = confint_qq$point[2,], y = object$dat),  mapping = aes(x=x,y=y), col = "grey") +
        #     geom_line(data = data.frame(x = confint_qq$overall[1,], y = object$dat),  mapping = aes(x=x,y=y), lty = 2) +
        #     geom_line(data = data.frame(x = confint_qq$overall[2,], y = object$dat),  mapping = aes(x=x,y=y), lty = 2) +
        #     geom_point()
        # }
      } else if(pl == "tmd"){
        xp <- qmod(p = txpos, scale = scale, shape = shape, family = object$family)
        yp <- dat
        xmap <- (xp + yp) / 2
        ymap <- yp - xp
        pl_list[["tmd"]] <-
          ggplot2::ggplot(data = data.frame(x = xmap, y = ymap),
                 mapping = ggplot2::aes(x = .data[["x"]],
                                         y = .data[["y"]])) +
          ggplot2::geom_hline(yintercept = 0, col = "gray") +
          ggplot2::geom_point() +
          ggplot2::labs(
               x = "average quantile",
               y = "quantile difference") +
          ggplot2::theme_classic()

      } else if(pl == "erp"){
        pl_list[["erp"]] <-
          ggplot2::ggplot(
            data = data.frame(
             y = ecdffun2(dat),
             x = ecdffun2(qmod(p = txpos,
                              scale = scale,
                              shape = shape,
                              family = object$family))),
          mapping = ggplot2::aes(
            x = .data[["x"]],
            y = .data[["y"]])) +
          ggplot2::geom_abline(intercept = 0, slope = 1, col = "gray") +
          ggplot2::geom_point() +
          ggplot2::labs(x = "theoretical quantiles",
               y = "empirical quantiles") +
          ggplot2::theme_classic()
      }
    }
    if(plot){
      lapply(pl_list, get("print.ggplot", envir = loadNamespace("ggplot2")))
    }
    return(invisible(pl_list))
  }
}



#' Uncertainty quantification for quantile-quantile plots
#' 
#' @param B number of bootstrap samples
#' @param dat vector of data
#' @param par parameter estimates of the model
#' @param lower lower bounds (truncation or lowest possible value)
#' @param upper upper bound for right-censoring or right-truncation
#' @param level level of the confidence intervals
#' @inheritParams nll_elife
#' @keywords internal
#' @export
#' @return a matrix with plotting positions for confidence intervals of quantile-quantile plots
uq1_qqplot_elife <-
  function(B = 9999L,
           dat,
           par,
           lower,
           upper,
           level = 0.95,
           type2 = c("none","ltrt","ltrc"),
           family = c("exp","gp","gomp","gompmake","weibull","extgp")
  ){
    n <- length(dat)
    family <- match.arg(family)
    if(missing(lower)){
      ltrunc <- 0
    }
    if(missing(upper)){
      rtrunc <- Inf
    }
    if(!missing(lower) && !missing(upper)){
      if(length(upper) != 1 && length(upper) != 1){
      stopifnot( "`upper` and `lower` should be vectors of the same length." = length(lower) == length(upper),
                 "`Length of data `dat` does not match vector of lower and upper bounds." = n == length(upper))
      }
    }
   stopifnot("`lower` should be smaller than `upper`." = isTRUE(all(lower <= upper)),
              "Number of bootstrap samples must be larger than what is prescribed by the level." = B >= 1/(1-level) - 1L)
    # parametric bootstrap samples
    # - simulate new data with the same sampling scheme
    # - estimate parameters of the distribution
    # - compute quantiles corresponding to plotting positions
    ppos <- matrix(NA, nrow = B, ncol = n)
    split_pars <- function(par, family){
    if(family == "gompmake"){
      scale <- as.numeric(par[-2])
      shape <- as.numeric(par[2])
    } else{
      scale <- as.numeric(par[1])
      shape <- as.numeric(par[-1])
    }
      return(list(scale = scale, shape = shape))
    }
    scale <- split_pars(par, family = family)$scale
    shape <- split_pars(par, family = family)$shape

    if(type2 == "none"){
      for(b in seq_len(B)){
        dat_boot <- samp_elife(n = n,
                               scale = scale,
                               shape = shape,
                               family = family,
                               type2 = "none")
        fit_boot <- fit_elife(time = dat_boot,
                              family = family)
        ppos[b,] <- qelife(p = ppoints(n),
                           scale = split_pars(fit_boot$par, family = family)$scale,
                           shape = split_pars(fit_boot$par, family = family)$shape,
                           family = family)

      }
    } else if(type2 == "ltrc"){
      for(b in seq_len(B)){
      boot_dat <- samp_elife(n = n,
                             scale = scale,
                             shape = shape,
                             lower = lower,
                             upper = upper,
                             family = family,
                             type2 = type2)
       fit_boot <- fit_elife(
        time = boot_dat$dat,
        thresh = 0,
        ltrunc = lower,
        event = !boot_dat$rcens,
        family = family,
        type = "right")
       np_boot <- np_elife(time = boot_dat$dat,
                event = !boot_dat$rcens,
                type = "right",
                ltrunc = lower)
       ppos[b,] <- qelife(p = n/(n+1)*np_boot$cdf(dat),
                          scale = split_pars(fit_boot$par, family = family)$scale,
                          shape = split_pars(fit_boot$par, family = family)$shape,
                          family = family)
      }
      } else if(type2 == "ltrt"){
        for(b in seq_len(B)){
          boot_dat <- samp_elife(n = n,
                                 scale = scale,
                                 shape = shape,
                                 lower = lower,
                                 upper = upper,
                                 family = family,
                                 type2 = type2)
          fit_boot <- fit_elife(time = boot_dat,
                                ltrunc = lower,
                                rtrunc = upper,
                                family = family)
          np_boot <- np_elife(time = boot_dat,
                            type = "interval",
                            event = rep(1L, n),
                            ltrunc = lower,
                            rtrunc = upper)$cdf
          xpos <- n / (n + 1) * (np_boot(dat) - np_boot(lower))/(np_boot(upper) - np_boot(lower))
          ppos[b,] <- qelife(p = pelife(q = lower,
                             scale = split_pars(fit_boot$par, family = family)$scale,
                             shape = split_pars(fit_boot$par, family = family)$shape,
                             family = family)*(1-xpos) +
            xpos * pelife(q = upper,
                         scale = split_pars(fit_boot$par, family = family)$scale,
                         shape = split_pars(fit_boot$par, family = family)$shape,
                         family = family),
            scale = split_pars(fit_boot$par, family = family)$scale,
            shape = split_pars(fit_boot$par, family = family)$shape,
            family = family)
        }
      }
    return(ppos)
    # return(boot::envelope(mat = ppos, level = level))
  }

# #' Approximate uncertainty for diagnostic plots
# #'
# #' Approximate the uncertainty by resampling parameters
# #' estimates using a normal approximation to the maximum
# #' likelihood estimators.
# #'
# #' @param B integer; number of simulations
# #' @param object an object of class \code{elife_par}
# #' @param logscale logical; if \code{TRUE} (default), compute the approximation using \eqn{\log(\sigma)} rather than \eqn{\sigma}
# #' @keywords internal
# uq2_qqplot_elife <- function(B = 1e4,
#                              object,
#                              logscale = TRUE){
#   stopifnot("Cannot compute approximation to sampling distribution of maximum likelihood estimators: missing `vcov` or `par` values" = !is.null(object$par) && !is.null(object$vcov))
#   par <- object$par
#   vcov <- object$vcov
#
# }

#' Plot empirical distribution function
#' 
#' @importFrom stats plot.ecdf
#' @export
#' @param x argument of class \code{elife_ecdf}
#' @return base R plot of the empirical distribution function
#' @param ... additional arguments for the plot
plot.elife_ecdf <- function(x, ...){
  args <- list(...)
  args$main <- ""
  args$x <- x
  args$ylab <- expression(F[n](x))
  do.call(plot.ecdf, args = args)
}
