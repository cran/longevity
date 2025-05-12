#' Turnbull's sets
#'
#' Given censoring and truncation set,
#' compute the regions of the real line that
#' get positive mass and over which the distribution
#' function is well-defined.
#'
#' @note The function adds the square root of the machine tolerance to left bounds of interval censored data
#' so they are open.
#'
#' @references Frydman, H. (1994). \emph{A Note on Nonparametric Estimation of the Distribution Function from Interval-Censored and Truncated Observations}, Journal of the Royal Statistical Society. Series B (Methodological) \bold{56}(1), 71-74.
#' @references Turnbull, B. W. (1976). \emph{The Empirical Distribution Function with Arbitrarily Grouped, Censored and Truncated Data.} Journal of the Royal Statistical Society, Series B \bold{38}, 290-295.
#' @export
#' @keywords internal
#' @return a matrix with two columns containing the \code{left} and \code{right} bounds Of Turnbull sets
#' @inheritParams npsurv
turnbull_intervals <- function(
    time,
    time2 = NULL,
    status,
    ltrunc = NULL,
    rtrunc = NULL
){
  stopifnot(length(time) > 0,
            length(status) == length(time),
            length(time) == length(time2))
  if(!isTRUE(all(is.finite(time2[status == 3L])))){
    stop("Invalid input: user specified interval censored observations, but \"time2\" is missing.")
  }
  lcens <- ifelse(status == 2L, -Inf, time)
  rcens <- ifelse(status %in% c(1L, 2L), time,
                  ifelse(status == 0L, Inf, time2))
  if(!is.null(rtrunc)){
    rtrunc <- as.vector(rtrunc[is.finite(rtrunc)])
    #don't put Inf, otherwise wrong...
  }
  if(length(rtrunc) == 0L){ # NULL has length zero
    rtrunc <- numeric(0)
  }
  if(!is.null(ltrunc)){
    ltrunc <- as.vector(ltrunc[is.finite(ltrunc)])
  }
  if(length(ltrunc) == 0L){
    ltrunc <- numeric(0)
  }
    # Frydman's correction with interval-censored and truncated
    # Without truncation, rtrunc and ltrunc are NULL
    result <- .turnbull_intervals(Lcens = lcens,
                                  Rcens = rcens,
                                  status = status,
                                  Ltrunc = ltrunc,
                                  Rtrunc = rtrunc)
  colnames(result) <- c("left", "right")
  return(result)
}



#' Nonparametric estimation of the survival function
#'
#' The survival function is obtained through the EM algorithm
#' described in Turnbull (1976); censoring and truncation are
#' assumed to be non-informative.
#' The survival function changes only at the \code{J} distinct
#' exceedances \eqn{y_i-u} and truncation points.
#'
#' The unknown parameters of the model are \eqn{p_j (j=1, \ldots, J)}
#' subject to the constraint that \eqn{\sum_{j=1}^J p_j=1}.
#' @inheritParams npsurv
#' @param thresh double thresh
#' @param tol double, relative tolerance for convergence of the EM algorithm
#' @param weights double, vector of weights for the observations
#' @param method string, one of \code{"em"} for expectation-maximization (EM) algorithm or \code{"sqp"} for sequential quadratic programming with augmented Lagrange multiplie method.
#' @param maxiter integer, maximum number of iterations for the EM algorithm
#' @param ... additional arguments, currently ignored
#' @return a list with elements
#' \itemize{
#' \item \code{cdf}: right-continuous \code{stepfun} object defined by probabilities
#' \item \code{time}: matrix of unique values for the Turnbull intervals defining equivalence classes; only those with non-zero probabilities are returned
#' \item \code{prob}: \code{J} vector of non-zero probabilities
#' \item \code{niter}: number of iterations
#' }
#' @useDynLib longevity, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @references Turnbull, B. W. (1976). \emph{The Empirical Distribution Function with Arbitrarily Grouped, Censored and Truncated Data.} Journal of the Royal Statistical Society. Series B (Methodological) 38(\bold{3}), 290–295.
#' @references Gentleman, R. and C. J. Geyer (1994). \emph{Maximum likelihood for interval censored data: Consistency and computation}, Biometrika, 81(\bold{3}), 618–623.
#' @references Frydman, H. (1994). \emph{A Note on Nonparametric Estimation of the Distribution Function from Interval-Censored and Truncated Observations}, Journal of the Royal Statistical Society. Series B (Methodological) \bold{56}(1), 71-74.
#' @export
#' @examples
#' set.seed(2021)
#' n <- 20L
#' # Create fake data
#' ltrunc <- pmax(0, runif(n, -0.5, 1))
#' rtrunc <- runif(n, 6, 10)
#' dat <- samp_elife(n = n,
#'                   scale = 1,
#'                   shape = -0.1,
#'                   lower = ltrunc,
#'                   upper = rtrunc,
#'                   family = "gp",
#'                   type2 = "ltrt")
#' npi <- np_elife(time = dat,
#'                 rtrunc = rtrunc,
#'                 ltrunc = ltrunc)
#' print(npi)
#' summary(npi)
#' plot(npi)
np_elife <- function(time,
                     time2 = NULL,
                     event = NULL,
                     type = c("right",
                              "left",
                              "interval",
                              "interval2"),
                     thresh = 0,
                     ltrunc = NULL,
                     rtrunc = NULL,
                     tol = 1e-12,
                     weights = NULL,
                     method = c("em", "sqp"),
                     arguments = NULL,
                     maxiter = 1e5L,
                     ...) {
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = np_elife, call = call, arguments = arguments)
    return(do.call(np_elife, args = arguments))
  }
  type <- match.arg(type)
  method <- match.arg(method)
  stopifnot(
    "Argument `thresh` should be positive." = min(thresh) >= 0,
    "Argument `thresh` should be a vector of length 1." = length(thresh) == 1L,
    "Argument `time` missing." = !missing(time),
    "Argument `time` should be a vector" = is.vector(time),
    "Argument `tol` should be of length one." = length(tol) == 1L,
    "Argument `tol` should be positive." = isTRUE(tol > 0),
    "Argument `tol` should be smaller." = isTRUE(tol < 1e-4)
  )
  survout <- .check_surv(time = time,
                         time2 = time2,
                         event = event,
                         type = type)
  time <- survout$time
  n <- length(time)
  time2 <- survout$time2
  status <- survout$status
  if(is.null(weights)){
    weights <- rep(1, n)
  } else{
    stopifnot("Weights must be positive and finite." = isTRUE(all(weights > 0)) & isTRUE(all(is.finite(weights))),
    "\"weights\" must be NULL or a vector of the same length as time."  = length(weights) == length(time)
    )
  }
  if(!is.null(ltrunc) && length(ltrunc) == 1L){
    ltrunc <- rep(ltrunc, length(time))
  }
  if(!is.null(rtrunc) && length(rtrunc) == 1L){
    rtrunc <- rep(rtrunc, length(time))
  }
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2,
                  FALSE,
                  time > thresh[1])
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    if(!is.null(time2)){
    time2 <- time2[ind] - thresh[1]
    }
    status <- status[ind]
    if(!is.null(ltrunc)){
      if(isTRUE(all(is.matrix(ltrunc),  ncol(ltrunc) > 1))){
        ltrunc <- apply(ltrunc[ind,] - thresh[1], 1:2, function(x){ pmax(0, x)})
      } else{
        ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
      }
    }
    if(!is.null(rtrunc)){
      if(isTRUE(all(is.matrix(rtrunc), ncol(rtrunc) > 1))){
        rtrunc <- apply(rtrunc[ind,] - thresh[1], 1:2, function(x){ pmax(0, x)})
      } else{
        rtrunc <- rtrunc[ind] - thresh[1]
      }
    }
  }
  # End of threshold shifting
  if(!is.null(ltrunc) & is.vector(ltrunc)){
    if(isTRUE(any(ltrunc > time, ltrunc > time2, na.rm = TRUE))){
      stop("Left-truncation must be lower than observation times.")
    }
  }
  if(!is.null(rtrunc) & is.vector(rtrunc)){
    if(isTRUE(any(rtrunc < time,
                  rtrunc < time2,
                  na.rm = TRUE))){
      stop("Right-truncation must be lower than observation times.")
    }
  }
  n <- length(time)
  cens <- !isTRUE(all(status == 1L,
                      na.rm = TRUE))
  if (!is.null(rtrunc)) {
    if(isTRUE(all(is.matrix(rtrunc), ncol(rtrunc) > 1))){
      stopifnot(
      "`rtrunc` should be the same length as `time`" =  nrow(rtrunc) == n
      )
    } else{
    stopifnot(
      "`rtrunc` should be a vector" = is.vector(rtrunc),
      "`rtrunc` should be the same length as `time`" =  length(rtrunc) == n
    )
    }
  }
  if (!is.null(ltrunc)) {
    if(isTRUE(all(is.matrix(ltrunc), ncol(ltrunc) > 1))){
      stopifnot(
        "`ltrunc` should be the same length as `time`" =  nrow(ltrunc) == n
      )
    } else{
      stopifnot(
      "`ltrunc` should be a vector" = is.vector(ltrunc),
      "`ltrunc` should be the same length as `dat`" =  length(ltrunc) == n
    )
    }
  }
  if(is.null(ltrunc) && is.null(rtrunc)){
    trunc <- FALSE
  } else{
    trunc <- TRUE
    stopifnot("`ltrunc` should be smaller than `rtrunc`" = isTRUE(all(ltrunc < rtrunc, na.rm = TRUE)))
  }

  if(length(event) == 1L){
    event <- rep(event, length(time))
  }
  if(is.matrix(ltrunc)){
    method <- "em"
  }
  # Dispatch methods
  if(method == "em"){
    return(
      npsurv(
        time = time,
        time2 = time2,
        event = event,
        type = type,
        ltrunc = ltrunc,
        rtrunc = rtrunc,
        weights = weights,
        status = status,
        thresh = thresh,
        arguments = NULL,
        maxiter = as.integer(maxiter))
      )
  # } else if(method == "emR"){
  # unex <- turnbull_intervals(
  #   time = time,
  #   time2 = time2,
  #   status = status,
  #   ltrunc = ltrunc,
  #   rtrunc = rtrunc)
  # # This treats interval censored times as semi-open, as in (L, R]
  #   time <- ifelse(status %in% c(0L, 3L), time + 1e-13, time)
  #   # time2 <- ifelse(status %in% c(2L, 3L), time2 - 1e-10, time2)
  # J <- nrow(unex)
  # if (J < 2) {
  #   stop("Only one interval: consider increasing the size of `delta`.")
  # }
  # if(J > 2500L){
  #   stop("The pure R implementation is currently limited to 2500 unique failure times for memory reasons.")
  # }
  # if(is.null(ltrunc)){
  #   ltrunc <- -Inf
  # }
  # if(is.null(rtrunc)){
  #   rtrunc <- Inf
  # }
  # # Check new definition of semi-infinite intervals are semi-closed
  # #time <- ifelse(is.finite(time2), time, time + 1e-8)
  # #time2 <- ifelse(is.finite(time), time2, time2 - 1e-8)
  # dummy1 <- matrix(NA, nrow = J, ncol = n)
  # dummy2 <- matrix(NA, nrow = J, ncol = n)
  # for(i in seq_len(J)){
  #   dummy1[i,] <- unex[i,1] >= time & unex[i,2] <= time2
  #   dummy2[i,] <- unex[i,1] >= ltrunc & unex[i,2] <= rtrunc
  # }
  # # Reduce the data to range
  # intervals1 <- apply(dummy1, 2,
  #                     function(x){
  #                       interv <- which(x)
  #                       if(length(interv) == 0L){
  #                         return(rep(J+1L, 2))
  #                       } else{
  #                         return(range(interv))
  #                       }})
  # intervals2 <- apply(dummy2, 2,
  #                    function(x){
  #                      interv <- which(x)
  #                      if(length(interv) == 0L){
  #                        return(rep(J+1L, 2))
  #                      } else{
  #                        return(range(interv))
  #                      }})
  # cens_lb <- intervals1[1,]
  # cens_ub <- intervals1[2,]
  #
  # trunc_lb <- intervals2[1,]
  # trunc_ub <- intervals2[2,]
  # p <- rep(1/J, J)
  # pnew <- rep(0, J)
  # C_mat <- D_mat <- matrix(0, nrow = n, ncol = J)
  # n_iter <- 0L
  # n_iter_max <- 1e4L
  # for (i in seq_len(n)) {
  #   # this doesn't need to be updated
  #   C_mat[i, cens_lb[i]] <- 1
  # }
  # while (n_iter < n_iter_max) {
  #   n_iter <- n_iter + 1L
  #   # E-step
  #   # Update C_mat (only relevant for right-censored data)
  #   if (!cens){
  #     for (i in seq_len(n)) {
  #         ij <- cens_lb[i]:cens_ub[i]
  #         if(!isTRUE(all.equal(ij, J + 1L))){
  #          C_mat[i, ij] <- p[ij] / sum(p[ij])
  #         }
  #     }
  #   }
  #   if(!trunc){
  #     # Update D_mat
  #     for (i in seq_len(n)) {
  #       ij <- trunc_lb[i]:trunc_ub[i]
  #       if(!isTRUE(all.equal(ij, J + 1L))){
  #        D_mat[i, -ij] <- p[-ij] / sum(p[ij])
  #       }
  #     }
  #   }
  #
  #   # M-step: maximize the log-likelihood
  #   pnew <- colSums(C_mat + D_mat)
  #   pnew <- pnew / sum(pnew)
  #   pnew[pnew < 1e-14 * n] <- 0
  # if(max(abs(p - pnew)) < tol && n_iter_max != n_iter) {
  #   n_iter_max <- n_iter + 1L
  # }
  #   p <- pnew
  # }
  # if(vcov){
  # # Compute the hessian matrix of the log-likelihood
  # uij <- C_mat + D_mat
  # Omega_mat <- Eta_mat <- matrix(FALSE, nrow = n, ncol = J)
  # for(i in 1:n){
  #   Omega_mat[i, cens_lb[i]:cens_ub[i]] <- TRUE
  #   Eta_mat[i, trunc_lb[i]:trunc_ub[i]] <- TRUE
  # }
  # sum_Omega_i_sq <- as.numeric(Omega_mat %*% p)^2
  # sum_Eta_i_sq <- as.numeric(Eta_mat %*% p)^2
  # infomat <- matrix(0, J-1, J-1)
  # for(j in 1:(J-1)){
  #   for(k in j:(J-1)){
  #   infomat[j,k] <- -sum((Omega_mat[,k]-Omega_mat[,J])*(Omega_mat[,j]-Omega_mat[,J])/sum_Omega_i_sq) +
  #     sum((Eta_mat[,k]-Eta_mat[,J])*(Eta_mat[,j]-Eta_mat[,J])/sum_Eta_i_sq)
  #   infomat[k,j] <- infomat[j,k]
  #   }
  # }
  # covmat <- try(solve(-infomat), silent = TRUE)
  # if(is.character(covmat)){
  #   covmat <- NULL
  # }
  # } else{
  #   covmat <- NULL
  # }
  # survfun <- .wecdf(x = unex[,2],
  #                   w = p,
  #                   type = "surv")
  # std.err <- NULL
  # if(vcov && !is.null(covmat)){
  #   if(is.matrix(covmat) & nrow(covmat) >= 1){
  #   std.err <- vapply(1:nrow(covmat), function(j){
  #     sum(covmat[j:nrow(covmat),j:nrow(covmat)])
  #     }, numeric(1))
  #   }
  # }
  # invisible(
  #   list(
  #      survfun = survfun,
  #      surv = 1-cumsum(p)[-length(p)],
  #      stdsurv = std.err,
  #      xval = unex,
  #      prob = p,
  #      vcov = covmat,
  #      niter = n_iter))
  } else if(method == "sqp"){
    if(is.null(ltrunc) & is.null(rtrunc)){
      trunc <- FALSE
    }
    if(is.null(ltrunc)){
      ltrunc <- numeric(0)
    }
    if(is.null(rtrunc)){
      rtrunc  <- numeric(0)
    }
    # Create equivalence classes
    unex <- turnbull_intervals(
      time = time,
      time2 = time2,
      status = status,
      ltrunc = ltrunc,
      rtrunc = rtrunc)
    # Compute limits of Turnbull sets (ranges for alpha, beta's)
    limits <- .censTruncLimits(
      tsets = unex,
      lcens = time,
      rcens = time2,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      cens = as.logical(cens),
      trunc = as.logical(trunc)) + 1L
    # Number of intervals
    J <- nrow(unex)
    if(J > 1){
    coptim <- Rsolnp::solnp(
      pars = rep(1/J, J-1),
      fun = function(par, ...){
        np_nll(par = par,
               cens_lb = limits[,1],
               cens_ub = limits[,2],
               trunc_lb = limits[,3],
               trunc_ub = limits[,4],
               cens = cens,
               trunc = trunc,
               weights = weights)},
      ineqLB = rep(0, J),
      ineqUB = rep(1, J),
      ineqfun = function(par, ...){ c(par, 1-sum(par))},
      LB = rep(0, J-1),
      UB = rep(1, J-1),
      control = list(delta = 1e-8,
                     tol = 1e-12,
                     trace = FALSE))
    } else{
      coptim <- list(pars = 1, convergence = TRUE, niter = 0)
    }
    prob <- c(coptim$pars, 1-sum(coptim$pars))
    invisible(structure(
      list(
        cdf = .wecdf(x = thresh + unex[,2], w = prob),
        xval = thresh + unex,
        prob = prob,
        thresh = thresh,
        convergence = coptim$convergence == 0,
        niter = coptim$outer.iter),
      class = c("elife_npar", "npcdf")
    ))
  }
}

#' @export
print.elife_npar <- function(x, ...){
  # Compute restricted mean
  height <- 1-x$cdf(x$xval[,1]-1e-10)
  width <- diff(c(x$thresh, x$xval[,1]))
  rmean <- sum(height * width) + x$thresh
  quants <- stats::quantile(x$cdf, c(0.75, 0.5, 0.25))
  cat("Nonparametric maximum likelihood estimator\n\n")
  cat("Routine", ifelse(x$convergence, "converged", "did not converge"), "\n")
  cat("Number of equivalence classes:", nrow(x$xval),"\n")
  if(x$cdf(x$xval[nrow(x$xval),2]) == 1){
    cat("Mean: ", rmean,"\n")
  } else{
    cat("Restricted mean at upper bound", x$xval[nrow(x$xval),2], ":", rmean,"\n")
  }
  cat("Quartiles of the survival function:", quants)
}

#' @export
summary.elife_npar <- function(object, ...){
  height <- 1-object$cdf(object$xval[,1]-1e-10)
  width <- diff(c(object$thresh, object$xval[,1]))
  rmean <- sum(height * width) + object$thresh
  quant <- as.numeric(quantile(object$cdf, probs = c(0.25, 0.5, 0.75)))
  # cat(paste0("NPMLE of CDF: ", nrow(object$xval), " equivalence classes with summary\n"))
  c("Min." = min(object$xval),
    "1st Qu." = quant[1],
    "Median" = quant[2],
    "Mean" = rmean,
    "3rd Qu." = quant[3],
    "Max." = max(object$xval))
}

#' @export
plot.elife_npar <- function(x, ...){
  if(is.null(x$cdf)){
    stop("Object should have a \"cdf\" slot.")
  }
  args <- list(...)
  main <- ifelse(is.null(args$main), "", args$main)
  xlab <- ifelse(is.null(args$xlab), "x", args$xlab)
  ylab <- ifelse(is.null(args$ylab), "distribution function", args$ylab)
 plot(x$cdf, main = main, xlab = xlab, ylab = ylab, bty = "l")
}

#' Marginal log likelihood function of the nonparametric multinomial with censoring and truncation
#' @param par vector of \code{D-1} parameters
#' @param cens_lb index of interval in which death occurs
#' @param cens_ub index of interval in which death occurs (if death is observed), or else the largest interval.
#' @param trunc_lb vector of largest index for the lower truncation
#' @param trunc_ub vector of smallest index for the upper truncation
#' @return a scalar, the negative log likelihood value
#' @keywords internal
np_nll <- function(par,
                   cens_lb,
                   cens_ub,
                   trunc_lb,
                   trunc_ub,
                   cens,
                   trunc,
                   weights) {
  pf <- c(par, 1 - sum(par))
  if (isTRUE(any(c(pf > 1, pf < 0)))) {
    return(1e8)
  }
  llp <- 0
  for (i in seq_along(cens_lb)) {
    if(cens){
    llp <- llp + weights[i]*log(sum(pf[cens_lb[i]:cens_ub[i]]))
    }
    if(trunc){
      llp <- llp - weights[i]*log(sum(pf[trunc_lb[i]:trunc_ub[i]]))
    }
  }
  return(-llp)
}


#' Weighted empirical distribution function
#'
#' @param x vector of length \code{n} of values
#' @param w vector of weights of length \code{n}
#' @param type string, one of distribution function (\code{dist}) or survival function (\code{surv})
#' @author Adapted from spatstat (c) Adrian Baddeley and Rolf Turner
#' @keywords internal
.wecdf <- function(x, w, type = c("dist","surv")){
  # Adapted from spatstat (c) Adrian Baddeley and Rolf Turner
  type <- match.arg(type)
  stopifnot(length(x) == length(w))  #also returns error if x is multidimensional
  nbg <- is.na(x)
  x <- x[!nbg]
  w <- w[!nbg]
  n <- length(x)
  w <- w / sum(w)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  ox <- order(x)
  x <- x[ox]
  w <- w[ox]
  vals <- sort(unique(x))
  xmatch <- factor(match(x, vals), levels = seq_along(vals))
  wmatch <- tapply(w, xmatch, sum)
  wmatch[is.na(wmatch)] <- 0
  if(type == "dist"){
    rval <- approxfun(vals, cumsum(wmatch),
                      method = "constant",
                      yleft = 0,
                      yright = 1,
                      f = 0,
                      ties = "ordered")
  } else{
    rval <- approxfun(vals, c(1,1-cumsum(wmatch)[-length(wmatch)]),
                      method = "constant",
                      yleft = 1,
                      yright = 0,
                      f = 1,
                      ties = "ordered")
  }
  class(rval) <- c("elife_ecdf", "ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  invisible(rval)
}


#' Check survival output

#' @inheritParams npsurv
#' @return a list with transformed inputs or an error
#' @export
#' @keywords internal
.check_surv <- function(time,
                        time2 = NULL,
                        event = NULL,
                        type = c("right",
                                 "left",
                                 "interval",
                                 "interval2")){
  stopifnot("User must provide a time argument" = !missing(time))
   time <- as.numeric(time)
   n <- length(time)
   type <- match.arg(type)
   if(is.null(event)){
     event <- 1L
   }
   #stopifnot("Event vector is missing" = !is.null(event))
   event <- as.integer(event)
   if(length(event) == 1L){
     event <- rep(event, n)
   }
   stopifnot("Event should be a vector of integers between 0 and 3" = isTRUE(all(event %in% 0:3)))
   if(!is.null(time2) && type %in% c("interval","interval2")){
     time2 <- as.numeric(time2)
     stopifnot("Both `time` and `time2` should be numeric vectors of the same length." = length(time2) == length(time))
   } else if(!is.null(time2)){
     stop("`time2` should only be provided for `type=interval`.")
   }
   if(type == "right"){
     #warning("Event should be provided, else all events are assumed to be observed.")
     if(is.null(event)){
       time2 <- time
       status <- rep(1L, length(time))
     } else{
       time2 <- ifelse(as.logical(event), time, NA)
       status <- ifelse(as.logical(event), 1L, 0L)
     }
   } else if(type == "left"){
     time2 <- time
     time <- ifelse(as.logical(event), time, NA)
     status <- ifelse(as.logical(event), 1L, 2L)
   } else if(type == "interval"){
     time2 <- ifelse(event == 0L, NA,
                     ifelse(event %in% c(1L,2L), time, time2))
     time <- ifelse(event == 2L, NA, time)
     status <- event
   } else if (type == "interval2") {
     stopifnot("`time2` must be a numeric vector." = is.numeric(time2),
               "`time` and `time2` must be two vectors of the same length" = length(time2) == length(time))
     time <- ifelse(is.finite(time), time, NA)
     time2 <- ifelse(is.finite(time2), time2, NA)
     unknown <- (is.na(time) & is.na(time2))
     if(any(unknown)){
       stop("Invalid data; time is not resolved")
     }
   }
  status <- ifelse(is.na(time), 2L,
                   ifelse(is.na(time2), 0L,
                          ifelse(time == time2, 1L, 3L)))
  return(list(time = ifelse(is.na(time), -Inf, time),
              time2 = ifelse(is.na(time2), Inf, time2),
              status = status))
}

#' Nonparametric maximum likelihood estimation for arbitrary truncation
#'
#' The syntax is reminiscent of the \link[survival]{Surv} function, with
#' additional vectors for left-truncation and right-truncation.
#'
#' @note Contrary to the Kaplan-Meier estimator, the mass is placed in the interval
#' [\code{max(time), Inf}) so the resulting distribution function is not deficient.
#'
#' @param time excess time of the event of follow-up time, depending on the value of event
#' @param event status indicator, normally 0=alive, 1=dead. Other choices are \code{TRUE}/\code{FALSE} (\code{TRUE} for death).
#' For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
#' Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have experienced an event.
#' @param time2 ending excess time of the interval for interval censored data only.
#' @param type character string specifying the type of censoring. Possible values are "\code{right}", "\code{left}", "\code{interval}", "\code{interval2}".
#' @param ltrunc lower truncation limit, default to \code{NULL}
#' @param rtrunc upper truncation limit, default to \code{NULL}
#' @param weights vector of weights, default to \code{NULL} for equiweighted
#' @param arguments a named list specifying default arguments of the function that are common to all \code{elife} calls
#' @param ... additional arguments passed to the functions
#' @seealso \code{\link[survival]{Surv}}
#' @return a list with components
#' \itemize{
#' \item \code{xval}: unique ordered values of sets on which the distribution function is defined
#' \item \code{prob}: estimated probability of failure on intervals
#' \item \code{convergence}: logical; \code{TRUE} if the EM algorithm iterated until convergence
#' \item \code{niter}: logical; number of iterations for the EM algorithm
#' \item \code{cdf}: nonparametric maximum likelihood estimator of the distribution function
#' }
#' @export
#' @examples
#' # Toy example with interval censoring and right censoring
#' # Two observations: A1: [1,3], A2: 4
#' # Probability of 0.5
#'
#' test_simple2 <- npsurv(
#'   time = c(1,4),
#'   time2 = c(3,4),
#'   event = c(3,1),
#'   type = "interval"
#' )
npsurv <- function(time,
                   time2 = NULL,
                   event = NULL,
                   type = c("right",
                            "left",
                            "interval",
                            "interval2"),
                   ltrunc = NULL,
                   rtrunc = NULL,
                   weights = NULL,
                   arguments = NULL,
                   ...){
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = npsurv, call = call, arguments = arguments)
    return(do.call(npsurv, args = arguments))
  }
  stopifnot("`time` must be a numeric vector." = is.numeric(time))
  args <- list(...)
  thresh <- ifelse(is.null(args$thresh), 0, args$thresh)
  if(!is.null(args$status)){
    status <- args$status
  } else{
  survout <- .check_surv(time = time,
              time2 = time2,
              event = event,
              type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  }
  n <- length(time)
  stopifnot(length(status) == n)
  # Modified Jan. 2023 - pass truncation limits as vectors, removing NAs
  # that would appear with doubly truncated results
  unex <- turnbull_intervals(time = time,
                             time2 = time2,
                             status = status,
                             ltrunc = na.omit(as.numeric(ltrunc)),
                             rtrunc = na.omit(as.numeric(rtrunc)))
  J <- nrow(unex)
  if(isTRUE(any(status != 1L))){
    cens <- TRUE
  } else{
    cens <- FALSE
  }
  if(is.null(ltrunc) & is.null(rtrunc)){
    trunc <- FALSE
  } else{
    trunc <- TRUE
  }
  if(is.null(ltrunc)){
    ltrunc <- rep(-Inf, ifelse(trunc, length(time), 1))
  }
  if(is.null(rtrunc)){
    rtrunc <- rep(Inf, ifelse(trunc, length(time), 1))
  }
   if(!is.null(args$tol)){
    tol <- args$tol
    if(isTRUE(any(length(tol) != 1L,
                  tol > 1/J,
                  tol <= 0))){
      stop("Invalid argument \"tol\".")
    }
  } else{
    tol <- 1e-8
  }
  if(!is.null(args$zerotol)){
    zerotol <- args$zerotol
    if(isTRUE(any(!is.numeric(zerotol),
                  length(zerotol) != 1,
                  zerotol <= 0,
                  zerotol < tol))){
      stop("Invalid argument \"zerotol\".")
    }
  } else{
    zerotol <- 1e-10
  }
  if(!is.null(args$maxiter)){
    maxiter <- as.integer(args$maxiter)
    if(isTRUE(any(maxiter <= 1L,
                  maxiter > 1e7,
                  !is.finite(maxiter)))){
      stop("Invalid argument \"maxiter\": must be a finite integer no bigger than 1e7L.")
    }
  } else{
    maxiter <- 1e5L
  }
  if(is.null(weights)){
    weights <- rep(1, n)
  } else{
    stopifnot("Weights must be positive and finite." = isTRUE(all(weights > 0)) & isTRUE(all(is.finite(weights))))
  }
  # Check if double truncation
  if(is.matrix(ltrunc)){
    # Check dimensions
    stopifnot(is.matrix(ltrunc),
              ncol(ltrunc) == 2L,
              ncol(rtrunc) == 2L,
              nrow(ltrunc) == nrow(rtrunc),
              nrow(ltrunc) == n)
    if(!(isTRUE(all(time > ltrunc[,1] | time > ltrunc[,2])))){
      stop("Some observed times are not part of the observation windows.")
    }
    if(!(isTRUE(all(time2 < rtrunc[,1] | time2 < rtrunc[,2])))){
      stop("Some observed times are not part of the observation windows.")
    }
    if(!(isTRUE(all((time > ltrunc[,1] & time2 < rtrunc[,1]) |
                    (time > ltrunc[,2] & time2 < rtrunc[,2]))))){
      stop("Some censoring windows are not part of the observation windows.")
    }

    survfit_res <- .turnbull_em_dtrunc(
      tsets = unex,
      lcens = as.numeric(time),
      rcens = as.numeric(time2),
      ltrunc = as.matrix(ltrunc),
      rtrunc = as.matrix(rtrunc),
      cens = as.logical(cens),
      zerotol = as.numeric(zerotol),
      tol = as.numeric(tol),
      maxiter = as.integer(maxiter),
      trunc = as.logical(trunc),
      weights = as.numeric(weights)
    )
  } else{
  survfit_res <- .turnbull_em(
    tsets = unex,
    lcens = as.numeric(time),
    rcens = as.numeric(time2),
    ltrunc = as.numeric(ltrunc),
    rtrunc = as.numeric(rtrunc),
    cens = as.logical(cens),
    zerotol = as.numeric(zerotol),
    tol = as.numeric(tol),
    maxiter = as.integer(maxiter),
    trunc = as.logical(trunc),
    weights = as.numeric(weights)
  )
  }
  if(!survfit_res$conv){
    warning(paste("EM algorithm failed to converge after",
                  survfit_res$niter,
                  "iterations."))
  }
  probs <- as.numeric(survfit_res$p)
  invisible(structure(
    list(cdf = .wecdf(x = thresh + unex[,2], w = probs),
       xval = thresh + unex[probs > 0, ],
       prob = probs[probs > 0],
       thresh = thresh,
       niter = survfit_res$neval,
       convergence = survfit_res$conv,
       abstol = survfit_res$abstol),
    class = c("elife_npar", "npcdf")
  ))
}
