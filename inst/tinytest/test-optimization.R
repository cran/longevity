# # Optimization routines
# # Check that we get the same values for the MLE
# # as for "mev" package for GP and for special
# # cases of the exponential distribution where
# # the solution is available in closed-form
library(tinytest)

set.seed(1234)
n <- 1e5L
lower <- runif(n, max = 10)
upper <- runif(n, min = 10, max = 20)
# Generated interval truncated data
samp1 <- longevity::samp_elife(
  n = n,
  scale = 10,
  shape = 0.1,
  lower = lower,
  upper = upper,
  family = "gp",
  type2 = "ltrt")

fit1 <- longevity::fit_elife(
  time = samp1,
  ltrunc = lower,
  rtrunc = upper,
  family = "gp")

# Check that specifying truncation bounds that are 0/Inf
# does not impact the MLE returned
samp2 <- 2 + longevity::relife(n = 1e3, scale = 2, shape = -0.1, family = "gp")
fit_2a <- longevity::fit_elife(
  time = samp2,
  ltrunc = rep(0, length(samp2)),
  rtrunc = rep(Inf, length(samp2)),
  thresh = 2,
  family = "gp")
fit_2b <- longevity::fit_elife(
  time = samp2,
  thresh = 2,
  family = "gp")
expect_equivalent(fit_2a$par, fit_2b$par, tolerance = 1e-4)
if(requireNamespace("mev", quietly = TRUE)){
# Compare with MLE algorithm of Grimshaw
fit_2c <- mev::fit.gpd(
  xdat = samp2,
  threshold = 2)
expect_equivalent(fit_2a$par, fit_2c$par, tolerance = 1e-4)
}

# Check exponential data with known MLE
samp3 <- rexp(n = 100, rate = 0.5)
fit_3a <- longevity::fit_elife(
  time = samp3,
  family = "exp")
expect_equivalent(
  current = fit_3a$par,
   target = mean(samp3),
   info = "exponential")

samp4 <- longevity::samp_elife(
  n = 100,
  scale = 2,
  upper = 3,
  family = "exp",
  type2 = "ltrc")
fit4 <- longevity::fit_elife(
  time = samp4$dat,
  event = !samp4$rcens,
  type = "right",
  family = "exp")
expect_equal(
  as.numeric(fit4$par),
  sum(samp4$dat)/sum(!samp4$rcens),
  info = "exponential with right censoring")

