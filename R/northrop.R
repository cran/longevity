#' Piece-wise generalized Pareto distribution
#'
#' Density, distribution, quantile functions and random number generation from the
#' mixture model of Northrop and Coleman (2014), which consists of \code{m}
#' different generalized Pareto distributions over non-overlapping intervals
#' with \code{m} shape parameters and one scale parameter; the other scale parameters are
#' constrained so that the resulting distribution is continuous over the domain
#' and reduces to a generalized Pareto distribution if all of the shape parameters are equal.
#'
#' @name gppiece
#' @references Northrop & Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n sample size
#' @param scale positive value for the first scale parameter
#' @param shape vector of \code{m} shape parameters
#' @param thresh vector of \code{m} thresholds
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{\Pr[X \leq x]} otherwise, \eqn{\Pr[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, the values are returned on the log scale
#' @returns a vector of quantiles (\code{qgppiece}), probabilities (\code{pgppiece}), density (\code{dgppiece}) or random number generated from the model (\code{rgppiece})
NULL

#' @rdname gppiece
#' @export
#' @keywords internal
dgppiece <- function(x,
                     scale,
                     shape,
                     thresh,
                     log = FALSE
                     ){
  #parametrization is (sigma_1, xi_1, ..., xi_m)

  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
            )
  m <- length(shape)
  if(m == 1L){
    return(dgpd(x = x, loc = thresh, scale = scale, shape = shape, log = log))
  }
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  Ind <- as.integer(cut(x, c(thresh, Inf), include.lowest = TRUE))
  ldens <-  logp[Ind] - log(scale[Ind]) -
      ifelse(shape[Ind] == 0,
             (x-thresh[Ind])/scale[Ind],
      (1+1/shape[Ind])*log(pmax(0,1+shape[Ind]*(x-thresh[Ind])/scale[Ind]))
      )
   ldens[is.na(ldens)] <- -Inf
   ldens[is.na(x)] <- NA
   if(log){
     return(ldens)
   } else{
     return(exp(ldens))
   }
}

#' @rdname gppiece
#' @export
#' @keywords internal
pgppiece <- function(q,
                     scale,
                     shape,
                     thresh,
                     lower.tail = TRUE,
                     log.p = FALSE
){
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
  )
  m <- length(shape)
  if(m == 1L){
     return(pgpd(q = q, loc = thresh, scale = scale,
                 shape = shape, lower.tail = lower.tail, log.p = log.p))
  }
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  prob_interv <- -diff(c(exp(logp),0))
  Ind <- as.integer(cut(q, c(thresh, Inf), include.lowest = TRUE))
  cum_prob <- 1-exp(logp)
  F_upper <- c(vapply(seq_len(m-1), function(i){pgpd(q = w[i], scale = scale[i], shape = shape[i])}, numeric(1)), 1)
  ret_p <- cum_prob[Ind] +
    prob_interv[Ind]*(1 - ifelse(shape[Ind] == 0,
              exp((thresh[Ind] - q) / scale[Ind]),
              pmax(0, (1 + shape[Ind] * (q - thresh[Ind]) / scale[Ind]))^(-1/shape[Ind])
    )) / F_upper[Ind]
  # if(shape[m] < 0){
  #   # check that any point beyond the support gets cumulative value of 1
  #   outbound <- which(q[Ind == m] > thresh[m] - scale[m]/shape[m])
  #   if(length(outbound) > 0){
  #     ret_p[outbound] <- cum_prob
  #   }
  # }
  outzero <- which(q <= min(thresh))
  if(length(outzero) > 0){
    ret_p[outzero] <- 0
  }
  if(!lower.tail){
    ret_p <- 1-ret_p
  }
  if(log.p){
    return(log(ret_p))
  } else{
    return(ret_p)
  }
}

#' @rdname gppiece
#' @export
#' @keywords internal
qgppiece <- function(p,
                     scale,
                     shape,
                     thresh,
                     lower.tail = TRUE,
                     log.p = FALSE
){
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L,
            "\"log.p\" must be a logical." = (length(log.p) == 1L) & is.logical(log.p)
  )
  m <- length(shape)
  if(log.p){
    p <- exp(p)
  }
  if(m == 1L){
    return(qgpd(p = p, loc = thresh, scale = scale, shape = shape,
                lower.tail = lower.tail))
  }
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                )
            )
  cum_prob_interv <- c(1-exp(logp), 1)
  F_upper <- c(vapply(seq_len(m-1), function(i){
    pgpd(q = w[i], scale = scale[i], shape = shape[i])}, numeric(1)), 1)

  stopifnot("Some probabilities are smaller than zero or larger than one." = isTRUE(all(c(p >= 0, p <= 1))))
  if(!lower.tail){
    p <- 1 - p
  }
  Ind <- as.integer(cut(p, cum_prob_interv, include.lowest = TRUE))
  # Because of the division and numerical overflow, 1 gets mapped to *nearly* one
  ps <-  ifelse(p!=1, (p - (1-exp(logp)[Ind]))/(-diff(c(exp(logp),0))[Ind])*F_upper[Ind], 1)
  thresh[Ind] +
            ifelse(shape[Ind] == 0,
                  - scale[Ind] * log(1-ps),
                  scale[Ind] * ((1-ps)^(-shape[Ind]) - 1)/shape[Ind]
  )

}

#' @rdname gppiece
#' @export
#' @keywords internal
rgppiece <- function(n,
                     scale,
                     shape,
                     thresh
){
# Compute the probability p_j, then p_j - p_{j+1} is the
# probability of falling in the interval (v_j, v_{j+1})
# 1) simulate a multinomial vector
# 2) sample a truncated generalized Pareto on (v_j, v_{j+1})
# 3) add back the lowest threshold
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
  )
  m <- length(shape)
  if(m == 1L){
    return(thresh + relife(n = n, family = "gp", scale = scale, shape = shape))
  }
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  prob_interv <- -diff(c(exp(logp),0))
  #sample.int(m, prob = prob_interv, replace = TRUE, size = n)
  n_mixt <- rmultinom(n = 1, size = n, prob = prob_interv)
  samp <- vector(mode = "numeric", length = n)
  pos <- 1L
  for(j in seq_len(m-1)){
    if(n_mixt[j] > 0){
      samp[pos:(pos+n_mixt[j]-1)] <- thresh[j] +
        qgpd(runif(n_mixt[j])*pgpd(q = w[j], scale = scale[j], shape = shape[j]),
             scale = scale[j], shape = shape[j])
    }
    pos <- pos + n_mixt[j]
  }
  if(n_mixt[m] > 0){
    samp[pos:n] <- thresh[m] + qgpd(runif(n_mixt[m]), scale = scale[m], shape = shape[m])
  }
  return(samp[sample.int(n, n)])
}

#' Score test of Northrop and Coleman
#'
#' This function computes the score test
#' with the piecewise generalized Pareto distribution
#' under the null hypothesis that the generalized Pareto
#' with a single shape parameter is an adequate simplification.
#' The score test statistic is calculated using the observed information
#' matrix; both hessian and score vector are obtained through numerical differentiation.
#'
#' The score test is much faster and perhaps less fragile than the likelihood ratio test:
#' fitting the piece-wise generalized Pareto model is difficult due to the large number of
#' parameters and multimodal likelihood surface.
#'
#' The reference distribution is chi-square
#' @inheritParams fit_elife
#' @export
#' @param thresh a vector of thresholds
#' @param test string, either \code{"score"} for the score test or \code{"lrt"} for the likelihood ratio test.
#' @return a data frame with the following variables:
#' \itemize{
#' \item \code{thresh}: threshold for the generalized Pareto distribution
#' \item \code{nexc}: number of exceedances
#' \item \code{score}: score statistic
#' \item \code{df}: degrees of freedom
#' \item \code{pval}: the p-value obtained from the asymptotic chi-square approximation.
#' }
#' @examples
#' \donttest{
#' set.seed(1234)
#' n <- 100L
#' x <- samp_elife(n = n,
#'                 scale = 2,
#'                 shape = -0.2,
#'                 lower = low <- runif(n),
#'                 upper = upp <- runif(n, min = 3, max = 20),
#'                 type2 = "ltrt",
#'                 family = "gp")
#' test <- nc_test(
#'   time = x,
#'   ltrunc = low,
#'   rtrunc = upp,
#'   thresh = quantile(x, seq(0, 0.5, by = 0.1)))
#' print(test)
#' plot(test)
#' }
nc_test <- function(time,
                    time2 = NULL,
                    event = NULL,
                    thresh = 0,
                    ltrunc = NULL,
                    rtrunc = NULL,
                    type = c("right", "left", "interval", "interval2"),
                    weights = rep(1, length(time)),
                    test = c("score", "lrt"),
                    arguments = NULL,
                    ...){
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = nc_test, call = call, arguments = arguments)
    return(do.call(nc_test, args = arguments))
  }

  # Exclude doubly interval truncated data
  if(isTRUE(all(is.matrix(ltrunc),
                is.matrix(rtrunc),
                ncol(ltrunc) == ncol(rtrunc),
                ncol(rtrunc) == 2L))){
    stop("Doubly interval truncated data not supported")
  }
  test <- match.arg(test)
  stopifnot("Threshold is missing" = !missing(thresh))
  nt <- length(thresh) - 1L
  thresh <- sort(unique(thresh))
  stopifnot("Threshold should be at least length two" = nt >= 1L)
  res <- data.frame(thresh = numeric(nt),
                    nexc = integer(nt),
                    stat = numeric(nt),
                    df = integer(nt),
                    pval = numeric(nt))
  for(i in seq_len(nt)){
    fit0 <- fit_elife(time = time,
                      time2 = time2,
                      event = event,
                      thresh = thresh[i],
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      type = type,
                      family = "gp",
                      weights = weights,
                      export = FALSE)
    if(test == "score"){
      if (!requireNamespace("numDeriv", quietly = TRUE)) {
        stop(
          "Package \"numDeriv\" must be installed to use this function.",
          call. = FALSE
        )
      }
      score0 <- try(numDeriv::grad(func = function(x){
        nll_elife(par = x,
                  time = time,
                  time2 = time2,
                  event = event,
                  thresh = thresh[i:length(thresh)],
                  type = type,
                  ltrunc = ltrunc,
                  rtrunc = rtrunc,
                  family = "gppiece",
                  weights = weights
        )},
        x = c(fit0$par['scale'],
              rep(fit0$par['shape'], nt - i + 2L))
      ))
      hess0 <- try(numDeriv::hessian(func = function(x){
        nll_elife(par = x,
                  time = time,
                  time2 = time2,
                  event = event,
                  thresh = thresh[i:length(thresh)],
                  type = type,
                  ltrunc = ltrunc,
                  rtrunc = rtrunc,
                  family = "gppiece",
                  weights = weights
        )},
        x = c(fit0$par['scale'],
              rep(fit0$par['shape'], nt - i + 2L))
      ))
      score_stat <- try(as.numeric(score0 %*% solve(hess0) %*% score0))
      if(!inherits(score_stat, "try-error")){
        res$stat[i] <- score_stat
        res$df[i] <- as.integer(nt - i + 1L)
        res$pval[i] <- pchisq(q = score_stat, df = nt - i + 1L, lower.tail = FALSE)
      }
    } else{
      fit1 <- try(fit_elife(time = time,
                            time2 = time2,
                            event = event,
                            thresh = thresh[i:length(thresh)],
                            ltrunc = ltrunc,
                            rtrunc = rtrunc,
                            type = type,
                            family = "gppiece",
                            weights = weights,
                            export = FALSE))
      if(!inherits(fit1, "try-error") & isTRUE(fit1$convergence)){
        anova_tab <- anova(fit0, fit1)
        res$stat[i] <- anova_tab[2,3]
        res$df[i] <- as.integer(anova_tab[2,4])
        res$pval[i] <- anova_tab[2,5]
      }
    }
    res$nexc[i] <- fit0$nexc
  }
  res$thresh <- thresh[-length(thresh)]
  class(res) <- c("elife_northropcoleman", "data.frame")
  attr(res, "test") <- test
  return(res)
}

#' P-value plot
#'
#' The Northrop-Coleman tests for penultimate models are
#' comparing the piece-wise generalized Pareto distribution
#' to the generalized Pareto above the lower threshold.
#'
#' @export
#' @param x an object of class \code{elife_northropcoleman}
#' @param plot.type string indicating the type of plot
#' @param plot logical; should the routine print the graph if \code{plot.type} equals \code{"ggplot"}? Default to \code{TRUE}.
#' @return a base R or ggplot object with p-values for the Northrop-Coleman test against thresholds.
#' @param ... additional arguments for base \code{plot}
plot.elife_northropcoleman <- function(x,
                                       plot.type = c("base", "ggplot"),
                                       plot = TRUE,
                                       ...){
  plot.type <- match.arg(plot.type)
  stopifnot("Invalid object type" = inherits(x, what = "elife_northropcoleman"),
            "Not enough thresholds to warrant a plot" = nrow(x) >= 2L,
            "Not enough fit to produce a plot" = sum(is.finite(x$pval)) >= 2)
  args <- list(...)
  xlab <- args$xlab
  if(is.null(xlab) | !is.character(xlab)){
    xlab <- "threshold"
  }
  ylab <- args$ylab
  if(is.null(ylab) | !is.character(ylab)){
    ylab <- "p-value"
  }
  type <- args$type
  if(is.null(type) | !is.character(type)){
    type <- "b"
  }
  testname <- switch(attributes(x)$test,
                     "score" = "score test",
                     "lrt" = "likelihood ratio test")
  if(plot.type == "ggplot"){
    if(!requireNamespace("ggplot2", quietly = TRUE)){
      warning("You must install `ggplot2`: switching to base R plot.")
      plot.type = "base"
    } else{
      g1 <- ggplot2::ggplot(data = x,
                            mapping = ggplot2::aes(x = .data[["thresh"]],
                                                   y = .data[["pval"]])) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(x = xlab,
                      y = ylab,
                      caption = testname) +
        ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),
                                    limits = c(0,1),
                                    expand = c(0,0.1)) +
        ggplot2::theme_classic()
      if(plot){
        get("print.ggplot", envir = loadNamespace("ggplot2"))(g1)
      }
      return(invisible(g1))
    }

  }
  if(plot.type == "base"){
    plot(x = x$thresh,
         y = x$pval,
         bty = "l",
         xlab = xlab,
         ylab = ylab,
         type = type,
         ylim = c(0,1),
         yaxs = "i")
    mtext(testname, side = 3, adj = 1)
  }
}


#' @export
print.elife_northropcoleman <- function(x, ..., digits = NULL, quote = FALSE, right = TRUE,
                                        row.names = TRUE, max = NULL, eps = 1e-3){
  n <- length(row.names(x))
  if (length(x) == 0L) {
    cat(sprintf(ngettext(n, "data frame with 0 columns and %d row",
                         "data frame with 0 columns and %d rows"), n), "\n",
        sep = "")
  }  else if (n == 0L) {
    print.default(names(x), quote = FALSE)
    cat(gettext("<0 rows> (or 0-length row.names)\n"))
  }  else {
    if (is.null(max)){
      max <- getOption("max.print", 99999L)
    }
    if (!is.finite(max)){
      stop("invalid 'max' / getOption(\"max.print\"): ",
           max)
    }
    omit <- (n0 <- max%/%length(x)) < n
    m <- as.matrix(format.data.frame(
      if (omit){
        x[seq_len(n0), , drop = FALSE]
      } else { x
      }, digits = digits, na.encode = FALSE))
    m[, ncol(m)] <- format.pval(x$pval, digits = 3, eps = eps, na.form = "")
    if (!isTRUE(row.names))
      dimnames(m)[[1L]] <- if (isFALSE(row.names)){
        rep.int("", if (omit){
          n0
        } else {n})
      } else {
        row.names
      }
    print(m, ..., quote = quote, right = right, max = max)
    if (omit)
      cat(" [ reached 'max' / getOption(\"max.print\") -- omitted",
          n - n0, "rows ]\n")
  }
  invisible(x)
}

#' @inherit nc_test
#' @export
#' @keywords internal
#' @return the value of a call to \code{nc_test}
nc_score_test <- function(time,
                          time2 = NULL,
                          event = NULL,
                          thresh = 0,
                          ltrunc = NULL,
                          rtrunc = NULL,
                          type = c("right", "left", "interval", "interval2"),
                          weights = rep(1, length(time))){
  .Deprecated("nc_test", package = "longevity",
              msg = "`nc_score_test` is deprecated. \nUse `nc_test` with argument `test=\"score\"` instead.")

  # This is provided for backward compatibility for the time being
  nc_test(time = time,
          time2 = time2,
          event = event,
          thresh = thresh,
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          type = type,
          weights = weights,
          test = "score")
}
