# Check nesting structure

library(longevity)
library(tinytest)
# set.seed(123)
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
  shape = 1,
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
            event = samp$rcens == FALSE,
            family = "exp")
fit_gomp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == FALSE,
            family = "gomp")
fit_gp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == FALSE,
            family = "gp")
fit_gompmake <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == FALSE,
            family = "gompmake")
fit_extgp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == FALSE,
            restart = TRUE,
            check = TRUE,
            family = "extgp")
fit_weibull <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == FALSE,
            family = "weibull")
fit_weibull_alt <-
  fit_elife(time = samp$dat,
            ltrunc = lower,
            family = "weibull")
fit_extweibull <-
  fit_elife(time = samp$dat,
            ltrunc = lower,
            event = samp$rcens == FALSE,
            family = "extweibull")
fit_beard <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == FALSE,
  family = "beard")
fit_beardmake <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == FALSE,
  family = "beardmake")
fit_perks <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == FALSE,
  family = "perks")
fit_perksmake <- fit_elife(
  time = samp$dat,
  type = "right",
  ltrunc = lower,
  event = samp$rcens == FALSE,
  family = "perksmake")
devs <-
c(deviance(fit_exp),
  deviance(fit_gp),
  deviance(fit_gomp),
  deviance(fit_gompmake),
  deviance(fit_weibull),
  deviance(fit_extgp))
length(devs) == length(unique(devs))

devdiff <- c(
  deviance(fit_exp) - deviance(fit_gp), # 1
  deviance(fit_exp) -  deviance(fit_weibull), #2
  deviance(fit_gomp) - deviance(fit_extgp), #3
  deviance(fit_perksmake) - deviance(fit_beardmake), #4
  deviance(fit_perks) - deviance(fit_beard), # 5
  deviance(fit_exp) - deviance(fit_gomp), # 6
  deviance(fit_gp) - deviance(fit_extgp) # 7

)

# Regular testing
if(devdiff[1] >= 0){
expect_equal(
  anova(fit_exp, fit_gp)[2,5],
  pchisq(devdiff[1], df = 1, lower.tail = FALSE))
} else if(devdiff[1] < -1e-04){
 expect_error(anova(fit_exp, fit_gp))
}
if(devdiff[2] >= 0){
expect_equal(
  anova(fit_exp, fit_weibull)[2,5],
  pchisq(devdiff[2], df = 1, lower.tail = FALSE))
}  else if(devdiff[2] < -1e-04){
  expect_error(anova(fit_exp, fit_weibull))
}
if(devdiff[3] >= 0){
expect_equal(
  anova(fit_gomp, fit_extgp)[2,5],
  pchisq(devdiff[3],
         df = 1, lower.tail = FALSE))
} else if(devdiff[3] < -1e-04){
  expect_error(anova(fit_gomp, fit_extgp))
}
if(devdiff[4] >= 0){
expect_equal(
  anova(fit_perksmake, fit_beardmake)[2,5],
  pchisq(devdiff[4],
         df = 1, lower.tail = FALSE))
} else if(devdiff[4] < -1e-04){
  expect_error(anova(fit_perksmake, fit_beardmake))
}

if(devdiff[5] >= 0){
expect_equal(
  anova(fit_perks, fit_beard)[2,5],
  pchisq(devdiff[5], df = 1, lower.tail = FALSE))
}  else if(devdiff[5] < -1e-04){
  expect_error(anova(fit_perks, fit_beard))
}

# Nonregular, one parameter
if(devdiff[6] >= 0){
expect_equal(
  anova(fit_exp, fit_gomp)[2,5],
  0.5*(devdiff[6] < 1e-4) +
    0.5*pchisq(devdiff[6],
         df = 1, lower.tail = FALSE)
  )
} else if(devdiff[6] < -1e-04) {
  expect_error(anova(fit_exp, fit_gomp))
}
if(devdiff[7] >= 0){
expect_equal(
  anova(fit_gp, fit_extgp)[2,5],
  0.5*(devdiff[7] < 1e-4) +
    0.5*pchisq(devdiff[7], df = 1, lower.tail = FALSE)
)
} else if(devdiff[7] < -1e-04){
  expect_error(anova(fit_gp, fit_extgp))
}
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

