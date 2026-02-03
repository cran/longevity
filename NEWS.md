# Version 1.3

## Changes

- Fix evaluation environments when arguments is used (#4)

# Version 1.2.1

## Changes

- Changes due to update of `ggplot2` to 4.0.0 (S3 classes not exported anymore)


# Version 1.2

## Changes

- Function `gp_endpt_prof` renamed `endpoint.profile`
- New plot `endpoint.tstab` for threshold stability plot, with associated plot/autoplot methods.
- Bug fix: change incorrect events in `npsurv`

# Version 1.1

## Changes: 

- New S3 methods `AIC`/`BIC` for models
- Fix `npmle` threshold, expose number of iterations
- S3 methods `summary` and `print` for objects of class `npelife` now returns restricted mean plus the threshold. The `summary` method returns the vector with a custom implementation, rather than using the method for ecdf.
- `fit_elife` now returns a more informative error message if no data exceed the threshold.
- `fit_elife` correctly recover from failure to compute standard errors for `exp` model

# Version 1.0.0  (Release date 2023-11-12)

Initial version
