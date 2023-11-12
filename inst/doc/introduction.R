## ----setup, eval = TRUE, echo = FALSE-----------------------------------------
library(longevity)

## -----------------------------------------------------------------------------
thresh0 <- 36525
data(dutch, package = "longevity")
dutch1 <- subset(dutch, ndays > thresh0 & !is.na(ndays) & valid == "A")

## -----------------------------------------------------------------------------
args <- with(dutch1, list(
  time = ndays,  # time vector
  ltrunc = ltrunc, # left truncation bound
  rtrunc = rtrunc, # right truncation
  thresh = thresh0, # threshold (model only exceedances)
  family = "gp")) # choice of parametric model

## -----------------------------------------------------------------------------
tstab_c <- tstab(
  arguments = args,
  family = "gp", # parametric model, here generalized Pareto
  thresh = 102:108 * 365.25, # overwrites thresh
  method = "wald", # type of interval, Wald or profile-likelihood
  plot = FALSE) # by default, calls 'plot' routine
plot(tstab_c, 
     which.plot = "shape", 
     xlab = "threshold (age in days)")

## -----------------------------------------------------------------------------
(m1 <- fit_elife(arguments = args, 
          thresh = 105 * 365.25,
          family = "gp",
          export = TRUE))
m0 <- fit_elife(arguments = args, 
          thresh = 105 * 365.25,
          family = "exp")
anova(m1, m0)
plot(m1, which.plot = "qq")

