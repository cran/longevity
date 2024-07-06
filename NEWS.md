
# Version 1.1.0

## Changes: 

- S3 methods `summary` and `print` for objects of class `npelife` now returns restricted mean plus the threshold. The `summary` method returns the vector with a custom implementation, rather than using the method for ecdf.
- `fit_elife` now returns a more informative error message if no data exceed the threshold.
- `fit_elife` correctly recover from failure to compute standard errors for `exp` model

# Version 1.0.0  (Release date 2023-11-12)

Initial version
