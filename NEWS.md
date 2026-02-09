# fmrihrf 0.2.0

## New Features

* New `hrf_boxcar()` function for simple boxcar (step function) HRFs with optional normalization.
* New `hrf_weighted()` function for arbitrary weighted-window HRFs with constant or linear interpolation.
* `regressor()` now accepts a list of HRF objects for trial-varying HRF designs.
* New `plot.Reg()` method for visualizing regressor objects.
* New `plot_regressors()` for comparing multiple regressors on one plot (ggplot2 or base R).
* New `plot_hrfs()` for comparing multiple HRF shapes.
* New `print.HRF()` method for concise HRF summaries.

## Improvements

* Revised hemodynamic response and regressor vignettes.
* Expanded test suite for new HRF types and trial-varying regressors.

## Bug Fixes

* Fixed critical bug in `as_hrf()` where parameters stored in the `params` attribute were never used at evaluation time. The fix creates a closure that properly captures and applies parameters during evaluation.

# fmrihrf 0.1.0

* Initial CRAN release