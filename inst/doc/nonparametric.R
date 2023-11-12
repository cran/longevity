## -----------------------------------------------------------------------------
library(longevity)
left <- c(0,15,12,17,13,0,6,0,14,12,13,12,12,0,0,0,0,3,4,1,13,0,0,6,0,2,1,0,0,2,0)
right <- c(16, rep(Inf, 4), 24, Inf, 15, rep(Inf, 5), 18, 14, 17, 15,
           Inf, Inf, 11, 19, 6, 11, Inf, 6, 12, 17, 14, 25, 11, 14)
test <- np_elife(time = left,   # left bound for time
                 time2 = right, # right bound for time
                 type = "interval2", # data are interval censored
                 event = 3) # specify interval censoring, argument recycled

plot(test)

## -----------------------------------------------------------------------------
test$xval
print(test)

