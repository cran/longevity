# Check that simulated values are within bounds
# and fail with invalid parameter values

library(longevity)
library(tinytest)

# Check that simulated values are within bounds
n <- 1000L
# Generate random data from parametric families
family <- sample(eval(formals(relife)$family), size = 1)
args <- list(family = family)
# Check number of input parameters and sample scale, rate and shape
npar <- longevity:::.npar_elife(family = family, return_npar = TRUE)
if(npar[1] > 0){
  args$scale <- rexp(n = npar[1])
}
if(npar[2] > 0){
  args$rate <- rexp(n = npar[2])
}
if(npar[3] > 0){
  args$shape <- runif(n = npar[3])
}
# Create bounds on the distribution scale
plow1 <- runif(n, max = 0.1)
plow2 <- runif(n, min = 0.6, max = 0.8)
pupp2 <- runif(n, min = 0.8)
pupp1 <- runif(n, min = 0.3, max = 0.6)
# Compute corresponding quantiles
lower1 <- do.call(what = qelife, args = c(list(p = plow1, lower.tail = TRUE), args))
lower2 <- do.call(what = qelife, args = c(list(p = plow2, lower.tail = TRUE), args))
upper1 <- do.call(what = qelife, args = c(list(p = pupp1, lower.tail = TRUE), args))
upper2 <- do.call(what = qelife, args = c(list(p = pupp2, lower.tail = TRUE), args))
expect_true(isTRUE(all(lower1 < upper2)))
# First, simulate data from model with only single bounds
args$lower <- lower1
args$upper <- upper2
args$n <- n
s1 <- do.call(what = samp_elife, args = c(args, type2 = "ltrt"))
s2 <- do.call(what = samp_elife, args = c(args, type2 = "ltrc"))
expect_true(all(s1 > args$lower, s1 <= args$upper))
expect_true(all(s2$dat > args$lower,
                s2$dat[s2$rcens] == args$upper[s2$rcens],
                s2$dat[!s2$rcens] <  args$upper[!s2$rcens]))
# Similar checks for doubly interval truncated data
args$lower <- cbind(lower1, lower2)
args$upper <- cbind(upper1, upper2)
s3 <- do.call(what = samp_elife, args = c(args, type2 = "ditrunc"))
expect_true(all((s3 > lower1 & s3 < upper1) | (s3 > lower2 & s3 < upper2)))

# Check we get errors with invalid parameters
expect_error(relife(n = 10, scale = -1, shape = 0, family = "gp"))
expect_error(relife(n = 10, rate = rexp(2), shape = 0, family = "perksmake"))

