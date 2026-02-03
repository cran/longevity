# Log of the maximal data information prior
lprior_mdi_elife <- function(par,
                           family = c("exp",
                                      "gp",
                                      "gomp")){
  family <- match.arg(family)
  if(missing(par)){
    stop("Invalid parameter vector")
  }
  if(family == "exp"){
    stopifnot("Invalid \"par\": must be a scalar" = length(par) == 1L,
             "Parameter \"par\" must be finite"  = is.finite(par))
  } else {
      stopifnot("Invalid \"par\": must be a vector of length 2." = length(par) == 2L,
                "Parameter \"par\" must be finite"  = isTRUE(all(is.finite(par))))
    }
  obound <- switch(family,
                   exp = par[1] <= 0,
                   gp = par[1] <= 0 | par[2] < -1,
                   gomp = par[1] <=0 | par[2] <= 0)
  if(obound){
    return(-Inf)
  }
  if(family == "exp"){
    return(-log(par[1]))
  } else if(family == "gp"){
    return(-log(par[1]) - par[2] - 1)
  } else if(family == "gomp"){
    # stopifnot("Install package \"gsl\" to use \"lprior_mdi_elife\" with the Gompertz model.\n Try `install.packages(\"gsl\")`" = requireNamespace("gsl", quietly = TRUE))

    ## Define function internally to
    ## avoid an extra dependency
    expint_E1 <- function (x){

      polyval <- function (p, x)
      {
        if (length(x) == 0)
          return(c())
        if (length(p) == 0)
          return(0 * x)
        if (!is.vector(p, mode = "numeric") && !is.vector(p, mode = "complex"))
          stop("Argument 'p' must be a real or complex vector.")
        if (!is.vector(x) && !is.matrix(x))
          stop("Argument 'x' must be a real or complex matrix.")
        n <- length(p)
        y <- outer(x[1:length(x)], (n - 1):0, "^") %*% p
        dim(y) <- dim(x)
        return(y)
      }

      stopifnot(is.numeric(x) || is.complex(x))
      eps <- .Machine$double.eps
      x <- c(x)
      n <- length(x)
      y <- numeric(n)
      p <- c(-3.60269362633602e-09, -4.81953845214096e-07, -2.56949832211593e-05,
             -0.000697379085953419, -0.0101957352984579, -0.078118635592482,
             -0.301243289276271, -0.777380732573553, 8.26766195236648)
      polyv <- polyval(p, Re(x))
      k <- which(abs(Im(x)) <= polyv)
      if (length(k) != 0) {
        egamma <- 0.577215664901533
        xk <- x[k]
        yk <- -egamma - log(xk + (0+0i))
        j <- 1
        pterm <- xk
        term <- xk
        while (any(abs(term) > eps)) {
          yk <- yk + term
          j <- j + 1
          pterm <- -xk * pterm/j
          term <- pterm/j
        }
        y[k] <- yk
      }
      k <- which(abs(Im(x)) > polyv)
      if (length(k) != 0) {
        m <- 1
        xk <- x[k]
        nk <- length(xk)
        am2 <- numeric(nk)
        bm2 <- rep(1, nk)
        am1 <- rep(1, nk)
        bm1 <- xk
        f <- am1/bm1
        oldf <- rep(Inf, nk)
        j <- 2
        while (any(abs(f - oldf) > (100 * eps) * abs(f))) {
          alpha <- m - 1 + (j/2)
          a <- am1 + alpha * am2
          b <- bm1 + alpha * bm2
          am2 <- am1/b
          bm2 <- bm1/b
          am1 <- a/b
          bm1 <- 1
          f <- am1
          j <- j + 1
          alpha <- (j - 1)/2
          beta <- xk
          a <- beta * am1 + alpha * am2
          b <- beta * bm1 + alpha * bm2
          am2 <- am1/b
          bm2 <- bm1/b
          am1 <- a/b
          bm1 <- 1
          oldf <- f
          f <- am1
          j <- j + 1
        }
        y[k] <- exp(-xk) * f - (0+1i) * pi * ((Re(xk) < 0) &
                                                (Im(xk) == 0))
      }
      if (all(Im(y) == 0))
        y <- Re(y)
      return(y)
    }

    e1 <- exp(1/par[2])*expint_E1(x = 1/par[2])
    #The limit exp(1/b)*E1(1/b) = 0 as b -> 0
    return(-log(par[1]) + ifelse(is.na(e1),0,e1))
  }
}

#' Box-Cox transformation function
#'
#' Given a vector of parameters, apply the Box-Cox transformation.
#'
#' @export
#' @keywords internal
#' @param x vector of arguments
#' @param lambda vector of Box-Cox parameters
#' @return a vector of the same length as \code{x}
#' @importFrom "graphics" "par"
#' @importFrom "stats" "quantile"
boxcox_transfo <- function(x, lambda = rep(1, length(x))){
  stopifnot(length(x) == length(lambda),
            is.numeric(x),
            is.numeric(lambda))
  ifelse(lambda == 0, log(x), (x^lambda-1)/lambda)
}

# Posterior density for selected model using
# rust and the maximal data information prior

#' Log posterior distribution with MDI priors
#'
#' Log of the posterior distribution for excess lifetime
#' distribution with maximal data information priors.
#' @export
#' @return a vector proportional to the log posterior (the sum of the log likelihood and log prior)
#' @inheritParams nll_elife
lpost_elife <- function(par,
                        time,
                        time2 = NULL,
                        event = NULL,
                        type = c("right","left","interval","interval2"),
                        ltrunc = NULL,
                        rtrunc = NULL,
                        family = c("exp","gp","gomp"),
                        thresh = 0,
                        weights = rep(1, length(time)),
                        status = NULL,
                        arguments = NULL, ...){
  if(!is.null(arguments)){
    call <- match.call(expand.dots = FALSE)
    arguments <- check_arguments(func = lpost_elife, call = call, arguments = arguments)
    caller_env <- parent.frame()
    return(do.call(lpost_elife, args = arguments, envir = caller_env))
  }
  type <- match.arg(type)
  family <- match.arg(family)
  obound <- switch(family,
                   exp = par[1] <= 0,
                   gp = par[1] <= 0 | par[2] < -1,
                   gomp = par[1] <=0 | par[2] <= 0)
  if(obound){
    return(-1e20)
  }
loglik  <- -nll_elife(par = par,
          time = time,
          time2 = time2,
          event = event,
          type = type,
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          thresh = thresh,
          family = family,
          weights = weights,
          status = status)
  logprior <- lprior_mdi_elife(family = family, par = par)
  lpost <- loglik + logprior
  ifelse(is.finite(lpost), lpost, -1e20)
}

