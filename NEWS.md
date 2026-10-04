# fmrihrf 0.4.0

## Breaking changes

* `block_hrf(normalize = TRUE)` now computes each basis's absolute-peak
  normalization factor once, over the full blocked span, using the fixed
  reference-grid policy of `normalize_hrf(..., "unit_peak_per_basis")`.
  Previously each call normalized against its requested times: a scalar
  evaluation could become 1 even far from the peak, and the convolution and
  loop methods could disagree. Scalar, vector, and chunked queries now retain
  the same scale. Negative responses keep their sign and zero bases stay zero.
  Rebuild designs and refit analyses that used normalized blocked HRFs.

* Every positive `block_hrf()` width is now integrated, including widths below
  `precision`. Previously a short block returned the unintegrated impulse
  response, creating a discontinuity when width crossed the numerical step.
  `width = 0` retains the impulse convention. With `summate = FALSE`, the
  integral is divided by width; with finite `half_life`, decay attenuates the
  integrand but does not renormalize the averaging weights.

* `hrf_bspline()`, `hrf_bspline_generator()` and `HRF_BSPLINE` now return a
  basis that is zero at both ends of the span. Previously the last basis
  function equalled 1 at `span` and was cut to 0 beyond it, so any weight on
  it produced a step at the end of every response. The `N` functions are now
  the interior functions of a clamped basis with `N + 2` functions (the same
  construction as `fmrireg::estimate_hrf()`). Fitted coefficients and design
  matrices built from these bases change; `N` must be at least `degree - 1`.
  The tent basis (`hrf_tent_generator()`, `getHRF("tent")`, and
  `hrf_bspline(degree = 1)`) changes the same way: its tents now peak at
  interior knots and the basis is zero at both ends (with `N = 5` over 24 s,
  peaks at 4, 8, ..., 20 s instead of 4.8, ..., 24 s).

* `hrf_weighted(method = "constant")` now uses every weight: `n` weights fill
  `n` bins. With `width`, the window is split into `n` equal bins; with
  `times`, weight `i` covers `[times[i], times[i + 1])` and the last bin is as
  wide as the one before it. Previously the last weight only closed the
  window (and was returned at the single time `times[n]`). `normalize = TRUE`
  now makes all `n` weights sum to 1. The `method = "linear"` form is
  unchanged. Weight vectors that ended with a 0 as an end marker give the same
  values as before, with `span` one bin longer.

* Corrected the documentation of `normalize` in `hrf_boxcar()` and
  `hrf_weighted()`. A unit-area boxcar makes the GLM coefficient the
  integrated signal in the window, not its mean, and a weighted HRF's
  coefficient is a least-squares amplitude, not a weighted mean. Behaviour is
  unchanged.

## Plotting

* New `hrf_palette()` and matching ggplot2 scales `scale_colour_hrf()`,
  `scale_color_hrf()` and `scale_fill_hrf()`. The categorical palette has at
  least 3:1 contrast on white and near-black backgrounds and stays
  distinguishable under common colour-vision deficiencies; the ordered palette
  is for basis functions and parameter sweeps. All package plots use them.

* `plot_hrfs()` and `plot_regressors()` now plot every column of a basis set
  (`basis = "first"` restores the old behaviour), accept `layout = "stack"` for
  one panel per curve, and choose the palette with `palette`.
  `plot_hrfs()` gains `reference` for a dashed grey comparison curve.
  `plot_regressors()` accepts a `regressor_set()`, draws events as bars whose
  width is the event duration, and shows scan-time values with `samples`.

* `plot_regressors()`, `plot.Reg()` and `plot.FeatureReg()` evaluate with a
  precision matched to the plotting grid by default (`precision = NULL`), so
  sharp HRF edges are drawn where they occur.

* Inside knitr documents, `plot_hrfs()` and `plot_regressors()` return their
  result visibly and let knitr print the plot, like a ggplot object. This lets
  document themes add dark-mode figure versions. Assigning the result inside a
  chunk no longer draws it; print it explicitly.

* The base-graphics `plot()` methods use the same palette, mark onsets with
  ticks on the time axis, and draw multi-basis regressors one panel per basis
  function (`layout = "overlay"` for the old style).

* All vignette figures were redrawn with these helpers; the reconstruction
  section of the advanced vignette now fits basis weights to target HRFs.

## Improvements

* Improved vignette plots for narrow screens: compact legends, shared-axis
  panels for gamma libraries and reconstruction, and visible light/dark
  series colors. Comparison helpers respect the active ggplot2 theme and
  accept `draw = FALSE` for customization or vignette auto-printing.

* Fixed clipped peak annotations and respected onset transparency in event
  and feature-regressor plots. Clarified the spline example's 24-second
  support without changing its evaluated values.

* Updated the website and vignettes to albersdown 2.1.0. Vignettes now use
  its self-contained output format with light, dark, and phone-sized figures,
  replacing the copied theme assets. Figure resolution is limited to keep
  the source package compact.

* Added `feature_regressor()` for continuously sampled features (for example
  RMS energy). Each sample is a zero-order-hold bin of width dt, with
  optional pre-convolution centering and scaling. This is the continuous
  analogue of an amplitude-modulated event train, not a list of trials.

* Added fixed-scale HRF normalization with `normalize_hrf()` and the
  `hrf_norm` argument to `gen_hrf()` and `getHRF()`. Modes include Nilearn/SPM
  reference-grid scaling, unit peak, unit integral, and independent per-basis
  unit peaks. The existing `normalise_hrf()` and `normalize = TRUE` interfaces
  retain their per-basis unit-peak behavior.

* Factored the duplicated single-basis / multi-basis branching in
  `block_hrf()`, `evaluate.HRF()`, and `normalise_hrf()` into three internal
  helpers (`.weighted_combine`, `.normalise_result`, `.get_peaks`). No
  user-visible behavior change.

## Convolution engine consolidation

The package carried four evaluation engines. `conv` provided the best
speed/accuracy tradeoff in development benchmarks, so there is now one compiled
convolution engine.

* `method = "conv"` remains the default and is the only compiled convolution
  engine. It avoids FFT zero-padding overhead and is substantially faster than
  `loop` on large designs; exact speedups depend on design size, basis count,
  and `precision`.

* `method = "fft"` and `method = "Rconv"` are deprecated. They now evaluate via
  `conv` and warn. The FFT engine never repaid its zero-padding, because the
  HRF is short relative to the sampled design, and it was the only method that
  could fail outright (on an internal FFT size above ~1e7). `Rconv` was an R
  reimplementation of `conv` that silently fell back to `loop` whenever the
  grid was irregular or durations varied. Both implementations have been
  removed.

* `method = "loop"` is retained as the reference implementation. It can be
  more accurate at a given `precision` because event onsets are evaluated
  directly, and it remains the automatic fallback for regressors built from a
  list of per-event HRFs.

* Event onsets and block edges are no longer snapped to the internal grid.
  Each is placed at its exact sub-bin position, substantially reducing
  off-grid error at the default `precision = 0.33`. Blocks are projected onto
  the linear hat basis, preserving trapezoid accuracy while keeping the exact
  edge placement.

* `method = "loop"` no longer truncates blocked events at `hrf_span`. It now
  extends to `hrf_span + duration`, removing a residual error that did not
  shrink as `precision` decreased. Its support calculation no longer assumes
  a multi-point regular grid, so one-point and irregular grids work as well.
  Blocks whose onset precedes the requested grid are retained when their
  duration carries the response into it.

## Bug Fixes

* Let the generated Rcpp wrapper translate native exceptions after C++ stack
  unwinding, removing the remaining direct `Rf_error()` call (issue #42).

* Clarified that `summate = FALSE` computes a temporal average, whose peak
  can decrease for longer blocks; `normalize = TRUE` controls unit-peak
  scaling (issue #49).

* Preserved parameter metadata in closed HRF constructors and decorators without
  incorrectly warning that the captured parameters would be ignored.

* Fixed loop-based block regressors integrating the kernel beyond its declared
  span, where the convolution engine already truncated it. The restored SPMG
  undershoot exposed this existing tail discrepancy. Support is now applied
  before block integration on both paths.

* **Breaking:** corrected the SPMG canonical positive coefficient from `0.0833`
  to `1/120` and used the exact undershoot coefficient `1/(6*15!)`.
  `HRF_SPMG1`, `HRF_SPMG2`, and `HRF_SPMG3` now have the SPM double-gamma
  shape (about 8.9% undershoot relative to peak, previously 0.6%). This changes
  both raw scale and shape; existing analyses should be refitted.

* **Breaking:** the third column of `HRF_SPMG3` is now a genuine response
  dispersion derivative, holding positive-component mean and mass fixed and
  using SPM's `(h(d) - h(d + 0.01))/0.01` sign/step convention. Previously it
  was a second time derivative. `deriv(HRF_SPMG3, ...)` follows the corrected
  basis. Temporal derivatives remain analytic. These are raw, unorthogonalized
  continuous kernels; exact sampled SPM/Nilearn design compatibility is not
  implied. Existing `normalize`, `normalise_hrf()`, `normalize_hrf()`, and
  `hrf_norm` scaling policies are unchanged. Coefficients and derivative-based
  amplitude summaries must respect the selected scaling and basis geometry.

* Fixed issue #50: `feature_regressor()` now rejects matrix and array inputs
  instead of silently flattening them into one long feature. Pass one feature
  column at a time.

Addresses the defects reported in issue #45.

* **Breaking:** epoch (`duration > 0`) regressors are no longer scaled by
  `1/precision`. `evaluate()` on a `Reg` object summed microtime samples of a
  unit-height boxcar without a step-size factor, so the amplitude of every
  epoch regressor grew as the `precision` argument shrank -- at the default
  `precision = 0.33` an epoch column was inflated roughly 3.2x relative to the
  integral it was meant to approximate. The compiled engines now apply the same
  trapezoid quadrature `evaluate.HRF()` uses, so a block response is
  `amplitude * integral h(t - onset - u) du` over the block and converges as
  `precision` decreases. Point events (`duration = 0`) are unaffected.
  Fitted betas from epoch designs will change scale. Because this release also
  corrects response shapes and event timing, model fits and t-statistics cannot
  be assumed unchanged; rebuild design matrices and refit existing analyses.

* **Breaking:** `summate = FALSE` now takes effect on every evaluation method.
  It was silently ignored by the `conv` (the default), `fft`, and `Rconv`
  engines, which never received the flag; only `method = "loop"` honoured it.

* **Breaking:** `evaluate.HRF()` selects its impulse and block branches on
  `duration` alone rather than on `duration < precision`. A numerical setting
  could previously decide which model was evaluated: at `precision = 0.2` a
  duration of 0.19 was treated as an impulse and 0.20 as a block, a five-fold
  amplitude jump. Durations smaller than `precision` are now integrated as
  blocks.

* The accepted evaluation methods (`conv`, `fft`, `Rconv`, `loop`) now
  agree to within quadrature error on block regressors; the deprecated names
  `fft` and `Rconv` route to `conv`.

* `hrf_sine()` and `hrf_fourier()` no longer error on a scalar `t`. `vapply()`
  dropped the result to a plain vector when `length(t) == 1`, so the support
  mask failed with "incorrect number of subscripts on matrix". This was
  reachable from public API via `lag_hrf()`, `block_hrf()`, and
  `getHRF("fourier")`.

* `HRF_BSPLINE` and `hrf_bspline_generator()` now place interior knots from
  `span`, fixed when the object is constructed. `splines::bs()` was called
  without `knots=` and fell back to quantiles of whatever `t` was supplied, so
  the basis depended on the evaluation grid and disagreed with `hrf_bspline()`.
  Evaluating on a single time point produced a degenerate basis.

* The internal Daguerre basis normalizes each column against a fixed reference
  grid rather than against the caller's `t`. Including negative lags inflated
  the divisor and shrank every returned value.

* `hrf_gaussian()`, `hrf_mexhat()`, `hrf_inv_logit()`, and `hrf_lwu()` return 0
  for `t < 0`, as the other kernels already did. The regressor path masked
  negative lags itself, but `block_hrf()` and `lag_hrf()` sample the shape
  function directly and mixed in pre-onset values -- up to 2% of peak for a
  blocked Mexican hat. Note that `hrf_mexhat()` is discontinuous at `t = 0` as
  a result, since its formula is non-zero there.

# fmrihrf 0.3.1

## New Features

* Added a package-owned command line interface with installed `fmrihrf` wrapper,
  `fmrihrf_cli()`, and `install_cli()`.

## Improvements

* Removed an unused suggested dependency and tightened build-ignore rules for
  local check artifacts.

# fmrihrf 0.3.0

## Improvements

* Consolidated derivative method Rd aliases into parent help pages, reducing documentation redundancy.
* Added explicit `importFrom(utils, tail)` to avoid R CMD check NOTEs.

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
