# fmrihrf 0.3.0

## Improvements

* Consolidated derivative method Rd aliases into parent help pages, reducing documentation redundancy.
* Added explicit `importFrom(utils, tail)` to avoid R CMD check NOTEs.

## Bug Fixes

* Guarded `is.symbol()` before `as.character()` in internal eco atlas extraction to prevent errors on non-symbol inputs.

# fmrihrf 0.2.1

## Bug Fixes

* Fixed `hrf_bspline()` support handling so values for `t > span` (and `t < 0`) are zeroed instead of wrapping to onset-like values.
* Fixed `block_hrf()` block integration to include quadrature step-size scaling, making amplitudes stable across `precision`.
* Fixed `hrf_sine()` and `hrf_fourier()` to clamp support to `[0, span]` and return zero outside the modeled window.
* Fixed `normalise_hrf()` to use fixed normalization constants computed on the HRF support, avoiding data-dependent scaling across evaluation grids.
* Fixed `evaluate.HRF()` block-duration summation to use the same weighted integration scheme as `block_hrf()`.
* Fixed `evaluate.Reg(normalize = TRUE)` to normalize regressor outputs consistently across evaluation methods, including single-trial regressors with different durations.
* Fixed `block_hrf(summate = FALSE)` to return normalized block integration (for both single- and multi-basis HRFs) instead of the legacy pointwise-maximum behavior.

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
