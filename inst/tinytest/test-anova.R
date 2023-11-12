# Check nesting structure

library(longevity)
library(tinytest)
set.seed(123)
n <- rpois(n = 1, lambda = 150)
lower <- ifelse(runif(n) < 0.8,
                0,
                qgomp(runif(n, max = 0.9),
                      scale = 1,
                      shape = 1.5))
upper <- lower + rexp(n, rate = 1/10)
samp <- longevity::samp_elife(
  n = n,
  scale = 2,
  shape = 1.5,
  type2 = "ltrc",
  lower = lower,
  upper = upper,
  family = "gomp")
expect_true(
  isTRUE(all(samp$dat - lower > 0)))
# Fit the model with all families
fit_exp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "exp")
fit_gomp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gomp")
fit_gp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gp")
fit_gompmake <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gompmake")
fit_extgp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            restart = TRUE,
            family = "extgp")
fit_weibull <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "weibull")
fit_weibull_alt <-
  fit_elife(time = samp$dat,
            ltrunc = lower,
            family = "weibull")
fit_beard <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == 1L,
  family = "beard")
fit_beardmake <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == 1L,
  family = "beardmake")
fit_perks <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == 1L,
  family = "perks")
fit_perksmake <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == 1L,
  family = "perksmake")
devs <-
c(deviance(fit_exp),
  deviance(fit_gp),
  deviance(fit_gomp),
  deviance(fit_gompmake),
  deviance(fit_weibull),
  deviance(fit_extgp))
length(devs) == length(unique(devs))

# Regular testing
expect_equal(
  anova(fit_exp, fit_gp)[2,5],
  pchisq(deviance(fit_exp) - deviance(fit_gp),
         df = 1, lower.tail = FALSE))
expect_equal(
  anova(fit_exp, fit_weibull)[2,5],
  pchisq(deviance(fit_exp) -
           deviance(fit_weibull),
       df = 1, lower.tail = FALSE))
expect_equal(
  anova(fit_gomp, fit_extgp)[2,5],
  pchisq(deviance(fit_gomp) - deviance(fit_extgp),
         df = 1, lower.tail = FALSE))
expect_equal(
  anova(fit_perksmake, fit_beardmake)[2,5],
  pchisq(deviance(fit_perksmake) - deviance(fit_beardmake),
         df = 1, lower.tail = FALSE))
expect_equal(
  anova(fit_perks, fit_beard)[2,5],
  pchisq(deviance(fit_perks) - deviance(fit_beard),
         df = 1, lower.tail = FALSE))
# Nonregular, one parameter
expect_equal(
  anova(fit_exp, fit_gomp)[2,5],
  0.5*pchisq(deviance(fit_exp) - deviance(fit_gomp),
         df = 1, lower.tail = FALSE)
  )
expect_equal(
  anova(fit_gp, fit_extgp)[2,5],
  0.5*pchisq(deviance(fit_gp) - deviance(fit_extgp),
             df = 1, lower.tail = FALSE)
)
# Non-identifiable parameters in nested models
expect_error(anova(fit_exp, fit_gompmake))
expect_error(anova(fit_exp, fit_perksmake))
expect_error(anova(fit_exp, fit_beardmake))

# Invalid comparisons - non-nested models
expect_error(anova(fit_weibull, fit_gp))
expect_error(anova(fit_weibull, fit_extgp))
expect_error(anova(fit_weibull, fit_gompmake))
expect_error(anova(fit_weibull, fit_gomp))
expect_error(anova(fit_gomp, fit_gp))
expect_error(anova(fit_gompmake, fit_gp))
expect_error(anova(fit_gompmake, fit_extgp))

# Different survival configuration
expect_error(anova(fit_exp, fit_weibull_alt))

# Same distributions
expect_error(anova(fit_weibull, fit_weibull_alt))

