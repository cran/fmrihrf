#' Construct a Feature Regressor from a Sampled Time Series
#'
#' Creates a regressor from a continuously sampled feature (for example RMS
#' energy of an acoustic stimulus) by treating each sample as a zero-order-hold
#' bin of width \eqn{\Delta t}. The result is the Riemann-sum / ZOH approximation
#' of convolving the feature with an HRF, not a train of unit-mass impulses.
#'
#' Amplitude modulation on this sampling grid is the same linear model as
#' convolving the sampled series: if \eqn{x} is the (possibly centered) feature
#' and \eqn{H} is the convolution operator, the regressor is \eqn{Hx}. Zeros are
#' kept, so an all-TR / all-sample series stays an all-sample series. No
#' unmodulated companion regressor is added.
#'
#' Centering and scaling are applied to the feature **before** convolution.
#' For a whole-run series (`mask = NULL`),
#' \deqn{H(x-\mu\mathbf{1}) = Hx - \mu H\mathbf{1}.}
#' Away from the run edges, \eqn{H\mathbf{1}} is nearly constant (overlapping
#' HRFs sum to a plateau), so with a GLM intercept the centered and raw
#' columns are affinely equivalent. They differ by the HRF-length ramp of
#' \eqn{H\mathbf{1}} at the start and end of the run. Default `center = TRUE`
#' removes that boundary term; `center = FALSE` keeps it. `scale = "sd"`
#' only changes the feature's units (still pre-convolution). It is not
#' standardization of the final BOLD-space column after filtering.
#'
#' Use `mask` when the feature should be centered only during an on-period
#' (stimulus present, task on, and so on). Center and scale then use only the
#' on-samples, and off-mask samples stay 0:
#' \deqn{H[m(x-\mu_{\mathrm{on}})] = H(mx) - \mu_{\mathrm{on}} Hm.}
#' That is **not** an affine transform of the all-sample series. Pair it with
#' a separate boxcar for the on-period if you want presence and intensity as
#' two questions.
#'
#' Each sample is a zero-order-hold bin of width \eqn{\Delta t}, not a
#' unit-mass impulse. Evaluate with `precision` less than or equal to the
#' feature sampling interval. Compared with
#' `regressor(times, amplitude = values, duration = 0)`, the predicted BOLD
#' is smaller by about \eqn{\Delta t}.
#'
#' @md
#' @param values Numeric vector of feature samples. Matrix and array inputs are
#'   rejected; construct one feature regressor per column instead.
#' @param hrf The hemodynamic response function to convolve with the feature.
#'   Same types as [regressor()], except a list of per-event HRFs is not
#'   allowed. Defaults to `HRF_SPMG1`.
#' @param times Numeric vector of sample times in seconds, same length as
#'   `values`. Mutually exclusive with `dt`. Times must be strictly increasing
#'   and non-negative. The last bin width is the last inter-sample gap.
#' @param dt Positive sampling interval in seconds. Mutually exclusive with
#'   `times`. Sample times are `start + seq(0, by = dt, length.out = length(values))`.
#' @param start Start time in seconds used only when `dt` is supplied.
#'   Defaults to 0.
#' @param center Logical; if `TRUE` (default), subtract the mean of the
#'   (masked) samples before convolution. For an all-sample series this
#'   removes \eqn{\mu H\mathbf{1}}, including the HRF-length run-boundary ramp.
#' @param scale Character; `"none"` (default) leaves native units, `"sd"`
#'   divides by the standard deviation of the (masked) samples after centering.
#'   This z-scores the feature, not the convolved design column.
#' @param mask Optional logical vector the same length as `values`. Center and
#'   scale statistics are computed on `mask == TRUE` samples only; off-mask
#'   samples are set to 0. This is not equivalent to global all-sample
#'   centering.
#' @param span Temporal window in seconds for the HRF, passed to [regressor()].
#'   If `NULL`, the HRF's own span is used.
#'
#' @return An S3 object of class `c("FeatureReg", "Reg", "list")`. Evaluation
#'   uses the same convolution path as [regressor()].
#'
#' @examples
#' # 10 Hz envelope over 8 seconds, demeaned
#' dt <- 0.1
#' t <- seq(0, 8, by = dt)
#' rms <- abs(sin(2 * pi * t / 4))
#' feat <- feature_regressor(rms, dt = dt, hrf = HRF_SPMG1)
#'
#' grid <- seq(0, 12, by = 1)
#' y <- evaluate(feat, grid, precision = dt)
#'
#' # Same ZOH encoding as a duration-dt event regressor (no centering)
#' feat_raw <- feature_regressor(rms, dt = dt, center = FALSE, scale = "none")
#' ev <- regressor(t, HRF_SPMG1, duration = dt, amplitude = rms)
#' all.equal(evaluate(feat_raw, grid, precision = dt),
#'           evaluate(ev, grid, precision = dt))
#'
#' @seealso [regressor()], [evaluate()]
#' @export
feature_regressor <- function(values,
                              hrf = HRF_SPMG1,
                              times = NULL,
                              dt = NULL,
                              start = 0,
                              center = TRUE,
                              scale = c("none", "sd"),
                              mask = NULL,
                              span = NULL) {

  scale <- match.arg(scale)

  if (!is.null(dim(values))) {
    stop("`values` must be a numeric vector; matrix and array inputs are not ",
         "supported. Pass one feature column at a time.", call. = FALSE)
  }
  if (!is.numeric(values) || length(values) == 0L) {
    stop("`values` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("`values` must contain finite numeric values.", call. = FALSE)
  }

  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be a single logical value.", call. = FALSE)
  }

  has_times <- !is.null(times)
  has_dt <- !is.null(dt)
  if (has_times == has_dt) {
    stop("Supply exactly one of `times` or `dt`.", call. = FALSE)
  }

  n <- length(values)
  hrf_is_list <- is.list(hrf) && !inherits(hrf, "HRF")
  if (hrf_is_list) {
    stop("`hrf` cannot be a list of per-event HRFs for a feature regressor.",
         call. = FALSE)
  }

  if (has_dt) {
    if (!is.numeric(dt) || length(dt) != 1L || is.na(dt) || !is.finite(dt) ||
        dt <= 0) {
      stop("`dt` must be a single positive finite number.", call. = FALSE)
    }
    if (!is.numeric(start) || length(start) != 1L || is.na(start) ||
        !is.finite(start) || start < 0) {
      stop("`start` must be a single non-negative finite number.", call. = FALSE)
    }
    times <- start + seq(0, by = dt, length.out = n)
    duration <- rep(dt, n)
    dt_attr <- dt
  } else {
    times <- as.numeric(times)
    if (length(times) != n) {
      stop("`times` must have the same length as `values`.", call. = FALSE)
    }
    if (anyNA(times) || any(!is.finite(times))) {
      stop("`times` must contain finite numeric values.", call. = FALSE)
    }
    if (any(times < 0)) {
      stop("`times` must be non-negative.", call. = FALSE)
    }
    if (n < 2L) {
      stop("`times` must have at least two samples so the last bin width ",
           "can be inferred.", call. = FALSE)
    }
    if (is.unsorted(times, strictly = TRUE)) {
      stop("`times` must be strictly increasing.", call. = FALSE)
    }
    gaps <- diff(times)
    duration <- c(gaps, gaps[length(gaps)])
    dt_attr <- if ((max(gaps) - min(gaps)) <=
                   sqrt(.Machine$double.eps) * max(1, stats::median(gaps))) {
      gaps[[1L]]
    } else {
      NA_real_
    }
  }

  amplitude <- as.numeric(values)
  if (is.null(mask)) {
    stats_idx <- seq_len(n)
  } else {
    if (!is.logical(mask) || length(mask) != n || anyNA(mask)) {
      stop("`mask` must be a non-missing logical vector the same length as ",
           "`values`.", call. = FALSE)
    }
    stats_idx <- which(mask)
    if (length(stats_idx) == 0L) {
      stop("`mask` must contain at least one TRUE value.", call. = FALSE)
    }
  }

  if (center) {
    amplitude[stats_idx] <- amplitude[stats_idx] - mean(amplitude[stats_idx])
  }

  if (identical(scale, "sd")) {
    sdv <- stats::sd(amplitude[stats_idx])
    if (!is.finite(sdv) || sdv == 0) {
      stop("`scale = \"sd\"` requires a non-zero standard deviation after ",
           "centering.", call. = FALSE)
    }
    amplitude[stats_idx] <- amplitude[stats_idx] / sdv
  }

  if (!is.null(mask)) {
    amplitude[!mask] <- 0
  }

  span_arg <- if (is.null(span)) 40 else span
  out <- Reg(
    onsets = times,
    hrf = hrf,
    duration = duration,
    amplitude = amplitude,
    span = span_arg,
    summate = TRUE,
    drop_zero_amplitude = FALSE
  )

  class(out) <- c("FeatureReg", class(out))
  attr(out, "center") <- center
  attr(out, "scale") <- scale
  attr(out, "dt") <- dt_attr
  out
}


#' @export
#' @method print FeatureReg
print.FeatureReg <- function(x, ...) {
  n_samp <- length(x$onsets)
  hrf_name <- attr(x$hrf, "name") %||% "custom function"
  nb <- nbasis(x$hrf)
  dt_attr <- attr(x, "dt")
  centered <- isTRUE(attr(x, "center"))
  scale_attr <- attr(x, "scale") %||% "none"

  cat("-- FeatureReg ", paste(rep("-", 50), collapse = ""), "\n", sep = "")
  cat("   Samples:", n_samp, "\n")
  if (n_samp > 0L) {
    cat("   Time range:", round(min(x$onsets), 2), "s to",
        round(max(x$onsets), 2), "s\n")
    if (!is.null(dt_attr) && length(dt_attr) == 1L && !is.na(dt_attr)) {
      cat("   Sampling interval:", dt_attr, "s\n")
    } else {
      cat("   Sampling interval: irregular\n")
    }
  }
  cat("   Center:", centered, "\n")
  cat("   Scale:", scale_attr, "\n")
  if (n_samp > 0L) {
    cat("   Amplitude range:", round(min(x$amplitude), 2), "to",
        round(max(x$amplitude), 2), "\n")
  }
  cat("   HRF:", hrf_name, paste0("(", nb, " basis function",
                                  if (nb == 1L) ")" else "s)"), "\n")
  cat("   HRF span:", x$span, "s\n")

  invisible(x)
}


#' @rdname shift
#' @export
#' @method shift FeatureReg
#' @importFrom assertthat assert_that
shift.FeatureReg <- function(x, shift_amount, ...) {
  dots <- list(...)

  if (missing(shift_amount) && "offset" %in% names(dots)) {
    shift_amount <- dots$offset
  }

  if (missing(shift_amount)) {
    stop("Must supply `shift_amount` or `offset`.", call. = FALSE)
  }

  assert_that(is.numeric(shift_amount) && length(shift_amount) == 1,
              msg = "`shift_amount` must be a single numeric value")

  if (length(x$onsets) == 0L) {
    return(x)
  }

  shifted <- x$onsets + shift_amount
  if (any(shifted < 0)) {
    stop("Shifting would produce negative sample times.", call. = FALSE)
  }

  x$onsets <- shifted
  x
}


