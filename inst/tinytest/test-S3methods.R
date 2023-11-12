# Check that plots work with both base R and with ggplot2

library(tinytest)
library(longevity)

set.seed(1234)
n <- 100L
x <- samp_elife(n = n,
           scale = 2,
           shape = 0.5,
           lower = low <- runif(n),
           upper = upp <- runif(n, min = 3, max = 15),
           type2 = "ltrc",
           family = "weibull")
fit <- fit_elife(
  time = x$dat,
  ltrunc = low,
  event = !x$rcens,
  family = "weibull",
  export = TRUE)
fit2 <- fit_elife(
  time = x$dat,
  ltrunc = low,
  event = !x$rcens,
  family = "extweibull",
  export = TRUE)
anova(fit, fit2)
plot(fit, which.plot = "pp", plot.type = "base")
plot(fit, which.plot = "qq", plot.type = "base")
# Warning because erp == pp w/o censoring
plot(fit, which.plot = "erp", plot.type = "base")
plot(fit, which.plot = "exp", plot.type = "base")
plot(fit, which.plot = "tmd", plot.type = "base")

plot(fit, which.plot = "pp", plot.type = "ggplot")
plot(fit, which.plot = "qq", plot.type = "ggplot")
plot(fit, which.plot = "erp", plot.type = "ggplot")
plot(fit, which.plot = "exp", plot.type = "ggplot")
plot(fit, which.plot = "tmd", plot.type = "ggplot")

expect_true(isTRUE(all.equal(deviance(fit), as.numeric(logLik(fit))*-2)))
expect_true(length(x$dat) == nobs(fit))

set.seed(1234)
n <- 100L
x <- samp_elife(n = n,
                scale = 2,
                shape = 0.5,
                lower = low <- runif(n),
                upper = upp <- runif(n, min = 3, max = 15),
                type2 = "ltrt",
                family = "weibull")
fit <- fit_elife(
  time = x,
  ltrunc = low,
  event = upp,
  family = "weibull",
  export = TRUE)
plot(fit, which.plot = "pp", plot.type = "base")
plot(fit, which.plot = "qq", plot.type = "base")
# Warning because erp == pp w/o censoring
expect_warning(plot(fit, which.plot = "erp", plot.type = "base"))
plot(fit, which.plot = "exp", plot.type = "base")
plot(fit, which.plot = "tmd", plot.type = "base")

plot(fit, which.plot = "pp", plot.type = "ggplot")
plot(fit, which.plot = "qq", plot.type = "ggplot")
expect_warning(plot(fit, which.plot = "erp", plot.type = "ggplot"))
plot(fit, which.plot = "exp", plot.type = "ggplot")
plot(fit, which.plot = "tmd", plot.type = "ggplot")
