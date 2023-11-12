# #' Local excess lifetime model helper
# #'
# #' This function breaks down the likelihood contribution
# #' from the sample in intervals so that the total contribution
# #' of an observation is only counted once.
# #' @inheritParams hazard_elife
# #' @keywords internal
# # local_elife <- function(
# #     breakpoints,
# #     time,
# #     time2,
# #     event,
# #     type = c("right",
# #              "left",
# #              "interval",
# #              "interval2"),
# #     ltrunc = NULL,
# #     rtrunc = NULL,
# #     thresh = 0,
# #     arguments = NULL,
# #     ...){
# # if(!is.null(arguments)){
# #  call <- match.call(expand.dots = FALSE)
# #  arguments <- check_arguments(func = local_elife, call = call, arguments = arguments)
# #  return(do.call(local_elife, args = arguments))
# # }
# #   type <- match.arg(type)
# #   # Check entries
# #   # Check model
# #   # Check that the breakpoints
# #   #  lie within the range of the data
# #   stopifnot(breakpoints > 0)
# #   type <- match.arg(type)
# #   if(type == "ltrc"){
# #   # Lexis diagram: vertical slice
# #   ab <- as.logical(I(time > lower)*I(time < upper))
# #   #survive beyond age under consideration
# #   datu <- dat[ab]
# #   slowu <- pmax(lower, ltrunc[ab])
# #   rcensor <- rcens[ab]
# #   aboveupp <- datu > upper
# #   if(length(aboveupp) > 0){
# #     rcensor[aboveupp] <- TRUE
# #   }
# #   return(dat = datu, ltrunc = slowu, rcens = rcensor)
# #   }
# # }
#
#
# #' Profile likelihood for hazard of generalized Pareto
# #'
# #' This function returns the profile likelihood-based 95% confidence intervals
# #' of the hazard of the generalized Pareto distribution with left-truncated and right-censored observations.
# #' @inheritParams gpd_cens
# #' @return point estimate along with 95% confidence interval
# # prof_gpd_hazard_confint <- function(dat, rightcens, slow, thresh){
# #   confintmat <- matrix(0, nrow = length(thresh), ncol = 4)
# #   th <- thresh
# #   k <- 200L
# #   datd <- dat/365.25
# #   slowd <- slow/365.25
# #   b <- 0
# #   for(i in 1:length(thresh)){
# #     b <- b + 1L
# #     dth <- (th[i] - th[1])
# #     low <- (th[i] - th[1])
# #     upp <- ifelse(i > length(thresh)-1, 20, (th[i] - th[1])+1)
# #     par_est <- optim(par = c(0.7,-0.3), fn = function(x){gpd_intcens(par = c(1/x[1]-x[2]*dth, x[2]),
# #                                                                      dat = datd,
# #                                                                      rightcens = rightcens, low = low, upp = upp,
# #                                                                      slow = slowd, expo = FALSE)}, control = list(reltol = 1e-11))
# #     hazard_mle <- par_est$par[1]
# #     infoi <- numDeriv::hessian(x = par_est$par,
# #                                func =  function(x){gpd_intcens(par = c(1/x[1]-x[2]*dth, x[2]),
# #                                                                dat = datd,
# #                                                                rightcens = rightcens, low = low, upp = upp,
# #                                                                slow = slowd, expo = FALSE)})
# #     stderr.transfo <- sqrt(diag(solve(infoi)))[1]
# #     if(stderr.transfo < 1e-6){
# #       grid_psi <-  par_est$par[1] + seq(-0.4, 0.4, length = k)
# #     } else{
# #       grid_psi <-  par_est$par[1] + seq(-3*stderr.transfo, 4*stderr.transfo, length = k)
# #     }
# #     prof_vals <- xi_sigma_vals <- rep(0, k)
# #     for (j in 1:k) {
# #       opt_prof <- optimize(f = function(xi, haza, dat, slow, rightcens, dth){
# #         gpd_intcens(par = c(1/haza-xi*dth, xi), dat = dat, rightcens = rightcens,
# #                     slow = slow, low = low, upp = upp)},
# #         upper = min(1.1,1/(grid_psi[j]*dth)), lower = -1,
# #         haza = grid_psi[j], dat = datd,
# #         rightcens = rightcens, slow = slowd, dth = dth)
# #       xi_sigma_vals[j] <- opt_prof$minimum
# #       prof_vals[j] <- opt_prof$objective
# #     }
# #     infin <- which(prof_vals  > 1e10-1)
# #     if(length(infin) > 0){
# #       prof_vals[infin] <- NA
# #     }
# #     prof <- structure(list(psi = grid_psi, psi.max = hazard_mle[1],
# #                            pll = -prof_vals, maxpll = -par_est$value, std.err = stderr.transfo),
# #                       class = "eprof")
# #     confintmat[b,] <- c(confint_int(prof, parm = "profile",print = FALSE)[1:3], stderr.transfo)
# #   }
# #   return(confintmat)
# # }
#
