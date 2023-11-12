# Checking the nonparametric maximum
# likelihood estimator of Turnbull
# against the implementation.
# library(longevity)
# library(tinytest)

# [1] Toy example with interval censoring and right censoring
#
# Two observations: A1: [1,3], A2: 4
# Probability of 0.5
test_simple1 <- longevity::np_elife(
  time = c(1,4),
  time2 = c(3,4),
  event = c(3,1),
  type = "interval",
  alg = "sqp")
test_simple2 <- longevity::npsurv(
  time = c(1,4),
  time2 = c(3,4),
  event = c(3,1),
  type = "interval"
)

# interval::icfit(c(1,4), c(3,4))$pf

expect_equal(c(0.5,0.5), test_simple2$prob)
expect_equal(c(0.5,0.5), test_simple1$prob)
# Second example with right-censoring and a tie
# risk set at t=1 is 3, 1 failure
# risk set at t=2 is 1, 1 failure
test_simple3a <- longevity::npsurv(
  time = c(1, 1, 2),
  event = c(0, 1, 1),
  type = "right")
test_simple3b <- longevity::np_elife(
  time = c(1, 1, 2),
  event = c(0, 1, 1),
  type = "right",
  alg = "sqp")

# interval::icfit(Surv(time = c(1,1,2), event = c(0,1, 1))~1)$pf
expect_equal(c(1/3, 2/3),
                       test_simple3a$prob)
expect_equal(c(1/3, 2/3),
                       test_simple3b$prob)

# interval::icfit(L = c(1,1,2), R = c(1,Inf,2))

# This fails _problem with definition in Turnbull?_
set.seed(1234)
# [2] Left truncated data
n <- 100L
u <- rnorm(n = n, sd = 1)
y <- u + rexp(n = n, rate = 1/10)
# Bell-Lynden estimator
utimes <- sort(unique(y))
ratio <- sapply(utimes, function(qj){
  1 - sum(y == qj)/
    sum(u <= qj & qj <=y)
})
density_bl <- -diff(c(1, cumprod(ratio)))
## This is unexpectedly long
npmle_ltrunc_long <-
  longevity::npsurv(time = y,
         ltrunc = u,
         tol = 1e-12)
expect_equal(npmle_ltrunc_long$prob,
                       density_bl)

# [3] Left-truncated right-censored data
# Create fake data with ties
set.seed(1234)
n <- 100L
ltrunc <- pmax(0, rnorm(n, sd = 4))
# Since about half of lower truncation bounds are zero
# The following creates ties!
time <- rexp(10, rate = 1/2) + ltrunc
# 0 for right-censored, 1 for observed
rightcens <-
  sample(c(0L, 1L),
         size = n,
         replace = TRUE,
         prob = c(0.2, 0.8))
# But some tied times are censored!
#time <- time + ifelse(rightcens == 0L, 1e-10, 0)
# The ECDF is supported on unique failure times
tun <- sort(unique(time[rightcens == 1L]))
# Product limit estimator
wgt <- sapply(tun, function(ti) {
  1 - sum(time == ti & rightcens == 1L) /
    sum(ltrunc <= ti & ti <= time)
  # & ifelse(rightcens == 1L, TRUE, ti < time))
})
# Test the function from the 'longevity' package
npmle_ltruncrcens_long <-
  longevity::npsurv(
    time = time,
    event = rightcens,
    type = "right",
    ltrunc = ltrunc
  )
npmle_ltruncrcens_long2 <-
  longevity::np_elife(
    time = time,
    event = rightcens,
    type = "right",
    alg = "sqp",
    ltrunc = ltrunc
  )
# Tsai and Jewell product limit estimator
probs <- -diff(c(1, cumprod(wgt)))
# Extract intervals
xval <- npmle_ltruncrcens_long$xval
zeroprobinterv <- which(xval[,2]-xval[,1] > 0)
prob_em <- npmle_ltruncrcens_long$prob

if(length(zeroprobinterv) != 0L){
  xval <- xval[-zeroprobinterv,]
  prob_em <- prob_em[-zeroprobinterv]
}
# Expect only unique failure times
expect_true(
     isTRUE(all(xval[,1] == xval[,2])))
# Same right and left censoring
expect_equal(xval[,2],
                       tun)
# Same probabilities
expect_equal(prob_em,
                       probs)
if(requireNamespace("survival", quietly = TRUE)){

# [4] right censored data (Kaplan-Meier)
n <- rpois(n = 1, lambda = 100)
time <- rpois(n, lambda = 100)
# Create censoring
status <- ifelse(runif(n) < 0.5, 0, 1)
utime <- time[status == 0] # observed failure times
cstat <- longevity:::.check_surv(time,
                     event = status,
                     type = 'right')
unex <- longevity::turnbull_intervals(
  time = cstat$time,
  time2 = cstat$time2,
  status = cstat$status)
cLimits <- .censTruncLimits(
  tsets = unex,
  rcens = cstat$time2,
  lcens = cstat$time,
  ltrunc = c(rep(-Inf, length(cstat$time))),
  rtrunc = c(rep(Inf, length(cstat$time))),
  trunc = FALSE,
  cens = TRUE)
# Fit model
test_rcens_long <-
  longevity::npsurv(time = time,
         event = status,
         type = 'right')

# Compare to Kaplan-Meier
km <- survival::survfit(survival::Surv(time, status) ~ 1)
probsKM <- -diff(c(1, km$surv))

# Last observation is censored
lastCensored <- max(time) == max(utime)
if(!lastCensored){
  expect_equal(test_rcens_long$xval[,1],
                         test_rcens_long$xval[,2])
  expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[,1])
  expect_equal(test_rcens_long$prob,
                         probsKM[probsKM>0])
}  else{
  # largest observation is right-censored
  nm <- nrow(test_rcens_long$xval)
  expect_equal(test_rcens_long$xval[-nm, 1],
                         test_rcens_long$xval[-nm, 2])
  expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[-nm,1])
  expect_equal(test_rcens_long$prob[-nm],
                         probsKM[probsKM>0])
}

# Flip sign
status[time == max(time)] <- 1-status[time == max(time)]
test_rcens_long <-
  longevity::npsurv(time = time,
         event = status,
         type = 'right')
# Compare to Kaplan-Meier
km <- survival::survfit(survival::Surv(time, status) ~ 1)
probsKM <- -diff(c(1, km$surv))
if(lastCensored){
  expect_equal(test_rcens_long$xval[,1],
                         test_rcens_long$xval[,2])
  expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[,1])
  expect_equal(test_rcens_long$prob,
                         probsKM[probsKM>0])
} else {
  # largest observation is right-censored
  nm <- nrow(test_rcens_long$xval)
  expect_equal(test_rcens_long$xval[-nm, 1],
                         test_rcens_long$xval[-nm, 2])
  expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[-nm,1])
  expect_equal(test_rcens_long$prob,
                         c(probsKM[probsKM>0], 1-sum(probsKM)))
}
}
# [5] Example of Frydman (1994)
time <- c(2,4,6)
time2 <- c(3,5,8)
ltrunc <- c(1,1,2.5)

F_int <- longevity::turnbull_intervals(
  time = time,
  time2 = time2,
  ltrunc = ltrunc,
  status = rep(3,3))
test <- longevity::npsurv(
  time = time,
  time2 = time2,
  event = rep(3,3),
  type = "interval",
  ltrunc = ltrunc)
Fref <- cbind(time, c(ltrunc[3], time2[-1]))
expect_equal(F_int,
                       Fref,
                       check.attributes = FALSE)
expect_equal(test$prob,
                       c(0.5,0.25,0.25))

# [6] Interval censoring

# Icens only deals with finite bounds

if(requireNamespace("interval", quietly = TRUE)){
res1_icfit <- interval::icfit(L = c(1,1,4), R = c(3,2,Inf))$pf
res1_longev <- longevity::npsurv(time = c(1,1,4),
                  time2 = c(3,2,Inf),
                  type = "interval2")$prob
expect_equivalent(res1_icfit,res1_longev)

res2_icfit <- interval::icfit(L = c(1,1,2,3), R = c(2,2,3,4))$pf
res2_longev <- longevity::npsurv(time = c(1,2,3),
                  time2 = c(2,3,4),
                  type = "interval2",
                  weights = c(2,1,1))$prob
expect_equivalent(res2_icfit,res2_longev)

# AIDS example from Lindsey and Ryan

left <- c(0,15,12,17,13,0,6,0,14,12,13,12,12,0,0,0,0,3,4,1,13,0,0,6,0,2,1,0,0,2,0)
right <- c(16, rep(Inf, 4), 24, Inf, 15, rep(Inf, 5), 18, 14, 17, 15, Inf, Inf, 11, 19, 6, 11, Inf, 6, 12, 17, 14, 25, 11, 14)
test <- longevity::npsurv(time = left,
               time2 = right,
               type = "interval2")
# test2 <- interval::icfit(L = left, R = right, icfitControl = interval::icfitControl(maxit = 1e4, epsilon = 1e-15))
# expect_equivalent(test$prob, test2$pf, tolerance = 1e-5)

  test3 <- Icens::EM(cbind(left, right), maxiter = 1e4, tol = 1e-15)
  expect_equivalent(test$prob, test3$pf[test3$pf > 0], tolerance = 1e-7)

}
# [7] Example from Gentleman and Geyer (1994)
# Interval censoring
# Solution is 1/3 in (0,1], 1/3 in (1,2] and 1/3 in (2,3]
# whereas Turnbull may return 0.5,0,0.5
time = c(0,1,1,0,0,2)
time2 = c(1,3,3,2,2,3)
expect_equivalent(
  longevity::npsurv(time = time,
                  time2 = time2,
                  type = "interval2")$prob,
  rep(1/3,3))


if(requireNamespace("DTDA", quietly = TRUE)){
set.seed(2021)
n <- 200L
# Create fake data
ltrunc <- pmax(0, runif(n, -0.5, 1))
rtrunc <- runif(n, 6, 10)
dat <- samp_elife(n = n,
                  scale = 1,
                  shape = -0.1,
                  lower = ltrunc,
                  upper = rtrunc,
                  family = "gp",
                  type2 = "ltrt")
trunc_dtda <- DTDA::lynden(
  X = dat,
  U = ltrunc,
  V = rtrunc,
  boot = FALSE,
  error = 1e-15)
trunc_long <- longevity::np_elife(
  time = dat,
  ltrunc = ltrunc,
  rtrunc = rtrunc)

expect_equivalent(
  max(abs(trunc_long$prob - trunc_dtda$density)),
  0,
  tolerance = 1e-5)
}

