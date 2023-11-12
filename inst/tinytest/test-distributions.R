# Test functions in distributions.R

# Check parameter constraints
library(longevity)
# Create fake ordered data
ddata <- c(NA, -1, 0, sort(rexp(100)))
endpt_gp_neg <- c(4, 4 + 1e-8)
endpt_extgp_neg <- c(4 * log(2), 4 * log(2) + 1e-8)
# Evaluate distribution functions

models <- list(
  exp = list(dist = "exp",
               args = list(rate = 0.5)),
  weibull = list(dist = "weibull",
                 args = list(
                   scale = 2,
                   shape = 0.5)),
  gpdpos = list(dist = "gpd",
                     args = list(
                       scale = 2,
                      shape = 0.5)),
  gpdneg = list(dist = "gpd",
                     args = list(
                       scale = 2,
                      shape = -0.5)),
  gpd_exp = list(dist = "gpd",
                     args = list(
                       scale = 2,
                       shape = 0)),
  gomp = list(dist = "gomp",
                 args = list(
                   scale = 2,
                   shape = 0.5)),
  gomp_exp = list(dist = "gomp",
                  args = list(
                    scale = 2,
                    shape = 0)),
  gompmake = list(dist = "gompmake",
                  args = list(
                    scale = 2,
                    shape = 0.5,
                    lambda = 0.5)),
  gompmake_gomp = list(dist = "gompmake",
                        args = list(
                          scale = 2,
                          shape = 0.5,
                          lambda = 0)),
  gompmake_exp = list(dist = "gompmake",
                       args = list(
                         scale = 4,
                         shape = 0,
                         lambda = 0.25)),
  extgpneg = list(dist = "extgp",
                      args = list(
                        scale = 2,
                        shape1 = 0.5,
                        shape2 = -0.5)),
  extgpos = list(dist = "extgp",
                  args = list(
                    scale = 2,
                    shape1 = 0.5,
                    shape2 = 0.5)),
  extgp_gpdneg = list(dist = "extgp",
                  args = list(
                    scale = 2,
                    shape1 = 0,
                    shape2 = -0.5)),
  extgp_gomp = list(dist = "extgp",
                   args = list(
                     scale = 2,
                     shape1 = 0.5,
                     shape2 = 0)),
  extweibull_gpdpos = list(dist = "extweibull",
                    args = list(
                      scale = 2,
                      shape1 = 0.5,
                      shape2 = 1)),
  extweibullpos = list(dist = "extweibull",
                     args = list(
                       scale = 2,
                       shape1 = 0.25,
                       shape2 = 0.8)),
  extweibullneg = list(dist = "extweibull",
                       args = list(
                         scale = 2,
                         shape1 = -0.5,
                         shape2 = 0.5)),
  extweibull_exp = list(dist = "extweibull",
                        args = list(
                          scale = 2,
                          shape1 = 0,
                          shape2 = 1)),
  extweibull_weibull  = list(dist = "extweibull",
                             args = list(
                               scale = 2,
                               shape1 = 0,
                               shape2 = 0.5)),
  perks = list(dist = "perks",
               args = list(
                 rate = 0.5,
                 shape = 0.5
               )),
  perks_exp = list(dist = "perks",
               args = list(
                 rate = 0,
                 shape = 1
               )),
  beard = list(dist = "beard",
                   args = list(
                     rate = 0.5,
                     shape1 = 0.5,
                     shape2 = 1
                   )),
  beard_exp = list(dist = "beard",
               args = list(
                 rate = 0,
                 shape1 = 1,
                 shape2 = 1
               )),
  beard_gomp = list(dist = "beard",
                   args = list(
                     rate = 0.25,
                     shape1 = 0.5,
                     shape2 = 0
                   )),
  beard_perks = list(dist = "beard",
                    args = list(
                      rate = 0.5,
                      shape1 = 0.5,
                      shape2 = 1
                    )),
 perksmake = list(dist = "perksmake",
                  args = list(
                    rate = 0.5,
                    shape = 0.5,
                    lambda = 1
                  )),
 perksmake_perks = list(dist = "perksmake",
                  args = list(
                    rate = 0.5,
                    shape = 0.5,
                    lambda = 0
                  )),
 beardmake_perks = list(dist = "beardmake",
                  args = list(
                    rate = 0.5,
                    shape1 = 0.5,
                    shape2 = 1,
                    lambda = 0
                  )),
 beardmake_perksmake = list(dist = "beardmake",
                        args = list(
                          rate = 0.5,
                          shape1 = 0.5,
                          shape2 = 1,
                          lambda = 1
                        )),
 beardmake_beard = list(dist = "beardmake",
                       args = list(
                         rate = 0.5,
                         shape1 = 0.5,
                         shape2 = 1,
                         lambda = 0
                       )),
 beardmake_exp = list(dist = "beardmake",
                          args = list(
                           rate = 0,
                           shape1 = 0.25,
                           shape2 = 4,
                           lambda = 0.375
                          )),
 beardmake_gomp = list(dist = "beardmake",
                      args = list(
                        rate = 0.25,
                        shape1 = 0.5,
                        shape2 = 0,
                        lambda = 0
                      )),
 beardmake_gompmake = list(dist = "beardmake",
                      args = list(
                        rate = 0.25,
                        shape1 = 0.5,
                        shape2 = 0,
                        lambda = 0.5
                      )),
 gppieceneg = list(dist = "gppiece",
                   args = list(
                     scale = 0.2,
                     shape = c(2, -0.1, -0.3),
                     thresh = c(0,0.5,0.6)
                   )),
 gppiecepos = list(dist = "gppiece",
                   args = list(
                     scale = 0.2,
                     shape = c(-0.1, 0.2),
                     thresh = c(0,0.5)
                   )),
 gppiece_gpdpos = list(dist = "gppiece",
                   args = list(
                     scale = 2,
                     shape = 0.5,
                     thresh = 0
                   )))

cdf <- simplify2array(lapply(models, function(x){
  args <- x$args
  args$q <- ddata
  do.call(paste0("p", x$dist), args = args)
}))

logsurv <- simplify2array(lapply(models, function(x){
  args <- x$args
  args$q <- ddata
  args$lower.tail <- FALSE
  args$log.p <- TRUE
  do.call(paste0("p", x$dist), args = args)
}))


expect_equal(as.numeric(cdf[1,]),
             rep(NA_real_, length(models)),
            info = "Distribution function is NA for missing values.")
expect_equal(as.numeric(cdf[2:3,]),
             rep(0, 2*length(models)),
             info = "Distribution function is zero for negative values.")
# Distribution function is nondecreasing

nondecreasing <- apply(cdf[-1,], 2, diff)
expect_true(current = isTRUE(all(nondecreasing >=0)),
            info = "Distribution function is nondecreasing.")

# Distribution function maps to unit interval
expect_true(current = isTRUE(
    all(cdf[-1,] >= 0,
        cdf[-1,] <= 1)),
            info = "Distribution function maps to [0,1].")
# Check nested models
wnested <- which(grepl(names(models), pattern = "_"))
sub <- match(sub(x = names(models)[wnested], pattern = ".*_", replacement = ""), names(models))

# Check special cases
for(i in seq_along(wnested)){
  test <- expect_equivalent(
  as.numeric(cdf[-1,wnested[i]]),
  as.numeric(cdf[-1,sub[i]]),
             info = paste("CDF", i,  names(models)[wnested[i]]))
  if(!test){
    print(test)
  }
}

# Check log survival matches
expect_equal(1-exp(logsurv[-1,]),
             cdf[-1,],
            info = "log survival functions match.")

# Check that density is zero outside of the domain
dens <- simplify2array(lapply(models, function(x){
  args <- x$args
  args$x <- ddata
  do.call(paste0("d", x$dist), args = args)
}))

expect_equal(as.numeric(dens[2,]),
             rep(0, length(models)),
             info = "Density zero for neg values")
# Check density is NA for NA
expect_equal(as.numeric(dens[1,]),
             rep(NA_real_, length(models)))
# Check that density is non-negative
expect_true(isTRUE(all(dens[-1,] >= 0)),
            info = "Density non-negative")
# Density integrates to 1

for(i in seq_along(models)){
 mod <- function(x){
   dist <- models[[i]]$dist
   args <- models[[i]]$args
   args$x <- x
   do.call(paste0("d", dist), args = args)
 }
 integ <- integrate(mod,
                    lower = 0,
                    upper = 500,
                    subdivisions = 1e5)
 test <- expect_equal(integ$value, 1,
                      tolerance = 1e-3,
              info = paste("Density of",
                           models[[i]]$dist,
                           "(model", i ,")",
                           "integrates to 1."))
if(!test){
  print(test)
}
}

# Density is zero outside of support
expect_equal(dgpd(
  endpt_gp_neg,
  loc = 0,
  scale = 2,
  shape = -0.5
), c(0, 0),
info = "GP density at endpoint is zero with negative shape")

expect_equal(dextgp(
  endpt_extgp_neg,
  scale = 2,
  shape1 = 0.5,
  shape2 = -0.5
),
c(0, 0),
info = "extgp density at endpoint is zero with negative shape")


# Check that CDF is inverse of quantile function
unif <- ppoints(100)

for(i in seq_along(models)){
  dist <- models[[i]]$dist
  args <- models[[i]]$args
  args$p <- unif
  args$q <- do.call(paste0("q", dist), args = args)
  args$p <- NULL;
  p <- do.call(paste0("p", dist), args = args)
  test <- expect_equal(p, unif,
                       info = paste0("CDF is inverse of quantile function for ",
                                     models[[i]]$dist,
                                     " (model ", i ,")"),
                       tol = 1e-5)
  if(!test){
    print(test)
  }
}

# CHECK HAZARD FUNCTION
hmodels <- models[! names(models) %in% c("gppieceneg", "gppiecepos", "gppiece_gpdpos")]
wnested <- which(grepl(names(hmodels), pattern = "_"))
sub <- match(sub(x = names(hmodels)[wnested], pattern = ".*_", replacement = ""), names(hmodels))

# Check that hazard is zero outside of the domain
hazard <- simplify2array(lapply(hmodels, function(x){
  args <- x$args
  args$x <- ddata
  do.call(paste0("h", x$dist), args = args)
}))

# Check hazard is NA for NA or negative values
expect_equal(as.numeric(hazard[1:2,]),
             rep(NA_real_, time = 2*ncol(hazard)))
# Check that hazard is non-negative when available
expect_true(isTRUE(all(hazard[!is.na(hazard)] >= 0)),
            info = "Hazard non-negative")


for(i in seq_along(wnested)){
  test <- expect_equivalent(
    as.numeric(hazard[-1,wnested[i]]),
    as.numeric(hazard[-1,sub[i]]),
    info = paste("hazard", i,  names(models)[wnested[i]]))
  if(!test){
    print(test)
  }
}
