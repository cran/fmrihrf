#' Lag an HRF Object
#'
#' Creates a new HRF object by applying a temporal lag to an existing HRF object.
#'
#' @param hrf The HRF object (of class `HRF`) to lag.
#' @param lag The time lag in seconds to apply. Positive values shift the response later in time.
#'
#' @return A new HRF object representing the lagged function.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' lagged_spmg1 <- lag_hrf(HRF_SPMG1, 5)
#' # Evaluate at time 10; equivalent to HRF_SPMG1(10 - 5)
#' lagged_spmg1(10)
#' HRF_SPMG1(5)
lag_hrf <- function(hrf, lag) {
  assertthat::assert_that(inherits(hrf, "HRF"), msg = "Input 'hrf' must be an HRF object.")
  assertthat::assert_that(
    is.numeric(lag) && length(lag) == 1 && is.finite(lag),
    msg = "'lag' must be a single finite numeric value."
  )

  # Original attributes
  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")

  # Create the lagged function
  lagged_func <- function(t) {
    hrf(t - lag)
  }

  # Create new HRF object using as_hrf
  .as_closed_hrf(
    f = lagged_func,
    name = paste0(orig_name, "_lag(", lag, ")"),
    nbasis = orig_nbasis,
    span = orig_span + max(0, lag), # Increase span if lag is positive
    params = c(orig_params, list(.lag = lag)) # Add lag to params for bookkeeping
  )
}


#' Create a Blocked HRF Object
#'
#' Creates a new HRF object representing a response to a sustained (blocked)
#' stimulus by convolving the input HRF with a boxcar function of a given width.
#'
#' @param hrf The HRF object (of class `HRF`) to block.
#' @param width The width of the block in seconds. Every positive width is
#'   integrated, including widths smaller than `precision`. Zero retains the
#'   original impulse response.
#' @param precision The sampling precision in seconds used for the internal convolution (default: 0.1).
#' @param half_life The half-life of an optional exponential decay applied during the block (default: Inf, meaning no decay).
#' @param summate Logical; if TRUE (default), responses within the block are
#'   integrated (summed). If FALSE, the integrated response is divided by the
#'   block width to obtain a temporal average. With finite `half_life`, decay
#'   attenuates the integrand but does not change this divisor. This does not fix peak
#'   height: longer blocks can have smaller peaks. Use `normalize = TRUE` for
#'   unit-peak scaling.
#' @param normalize Logical; if TRUE, scale each basis independently by its
#'   peak absolute value on a fixed reference grid over the full blocked span
#'   (the original span plus `width`). Uses the same approximately 0.02-second
#'   reference spacing as [normalize_hrf()] with `"unit_peak_per_basis"`.
#'   The factors are computed once at construction, so scalar, vector, and
#'   chunked queries have the same scale. Signs are preserved and zero bases
#'   remain zero. Default: FALSE.
#'
#' @return A new HRF object representing the blocked function.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' blocked_spmg1 <- block_hrf(HRF_SPMG1, width = 5)
#' t_vals <- seq(0, 30, by = 0.5)
#' plot(t_vals, HRF_SPMG1(t_vals), type = 'l', col = "blue", ylab = "Response", xlab = "Time")
#' lines(t_vals, blocked_spmg1(t_vals), col = "red")
#' legend("topright", legend = c("Original", "Blocked (width=5)"), col = c("blue", "red"), lty = 1)
block_hrf <- function(hrf, width, precision = 0.1, half_life = Inf, summate = TRUE, normalize = FALSE) {
  assertthat::assert_that(inherits(hrf, "HRF"), msg = "Input 'hrf' must be an HRF object.")
  assertthat::assert_that(
    is.numeric(width) && length(width) == 1 && width >= 0 && is.finite(width),
    msg = "'width' must be a single non-negative finite numeric value."
  )
  assertthat::assert_that(
    is.numeric(precision) && length(precision) == 1 && precision > 0 && is.finite(precision),
    msg = "'precision' must be a single finite positive numeric value."
  )
  assertthat::assert_that(
    is.numeric(half_life) && length(half_life) == 1 && half_life > 0,
    msg = "'half_life' must be a single positive numeric value."
  )
  assertthat::assert_that(is.logical(summate) && length(summate) == 1, msg = "'summate' must be a single logical value.")
  assertthat::assert_that(is.logical(normalize) && length(normalize) == 1, msg = "'normalize' must be a single logical value.")

  # Original attributes
  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")

  # Create the blocked function
  blocked_func <- function(t) {
    if (width == 0) {
      # Exactly zero denotes an impulse; a positive width always denotes a
      # block, regardless of the numerical integration step.
      res <- hrf(t)
    } else {
      quad <- .block_offsets_weights(width, precision)
      hmat_list <- lapply(quad$offsets, function(offset) {
        decay_factor <- if (is.infinite(half_life)) 1 else exp(-log(2) * offset / half_life)
        hrf(t - offset) * decay_factor
      })
      res <- .weighted_combine(hmat_list, quad$weights,
                                nbasis = orig_nbasis, summate = summate)
    }

    return(res)
  }

  # Store parameters used for blocking
  block_params <- list(
      .width = width,
      .precision = precision,
      .half_life = half_life,
      .summate = summate,
      .normalize = normalize
  )
  
  # Create new HRF object using as_hrf
  blocked <- .as_closed_hrf(
    f = blocked_func,
    name = paste0(orig_name, "_block(w=", width, ")"),
    nbasis = orig_nbasis,
    span = orig_span + width, # Span increases by the block width
    params = c(orig_params, block_params) # Add block params for bookkeeping
  )
  if (normalize) {
    # Share the existing fixed-grid policy, but retain this constructor's
    # name and parameter metadata. Never scale against a caller's query grid.
    normalized <- normalize_hrf(blocked, "unit_peak_per_basis")
    attributes(normalized) <- attributes(blocked)
    return(normalized)
  }
  blocked
}


#' Normalize an HRF Object with a Fixed Scale
#'
#' Creates an HRF whose evaluations are divided by constants computed once on
#' a fixed reference grid. The resulting scale therefore does not depend on the
#' time points supplied in later calls.
#'
#' @param hrf An object of class `HRF`.
#' @param mode Normalization convention: `"spm"` divides by the sum of the
#'   canonical basis on the 1,600-point 0--32 second Nilearn reference grid;
#'   `"unit_peak"` divides every basis by the canonical basis peak;
#'   `"unit_integral"` divides every basis by the canonical basis trapezoidal
#'   integral; `"unit_peak_per_basis"` scales each basis independently; and
#'   `"none"` returns `hrf` unchanged.
#'
#' @return An `HRF` object with fixed normalization.
#' @details For multi-basis HRFs, `"spm"`, `"unit_peak"`, and
#'   `"unit_integral"` use one scalar computed from the first (canonical)
#'   column and apply it uniformly. This preserves the relative scale of
#'   derivative bases. Only `"unit_peak_per_basis"` rescales columns
#'   independently.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' spm_scaled <- normalize_hrf(HRF_SPMG1, "spm")
#' reference_grid <- seq(0, 32, length.out = 1600)
#' sum(spm_scaled(reference_grid))
#'
#' peak_scaled <- normalize_hrf(HRF_SPMG2, "unit_peak")
#' max(abs(peak_scaled(seq(0, 24, by = 0.01))[, 1]))
normalize_hrf <- function(
    hrf,
    mode) {
  assertthat::assert_that(
    inherits(hrf, "HRF"),
    msg = "Input 'hrf' must be an HRF object."
  )
  mode <- match.arg(mode, c(
    "spm", "unit_peak", "unit_integral", "unit_peak_per_basis", "none"
  ))
  if (identical(mode, "none")) {
    return(hrf)
  }

  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")
  param_names <- names(orig_params)
  decorator_params <- if (is.null(param_names)) {
    list()
  } else {
    orig_params[startsWith(param_names, ".")]
  }

  if (identical(mode, "spm")) {
    ref_grid <- seq(0, 32, length.out = 1600L)
  } else {
    ref_n <- max(as.integer(round(orig_span * 50)) + 1L, 2L)
    ref_grid <- seq(0, orig_span, length.out = ref_n)
  }
  ref_vals <- hrf(ref_grid)
  canonical <- if (is.matrix(ref_vals)) ref_vals[, 1L] else as.numeric(ref_vals)

  factor <- switch(
    mode,
    spm = sum(canonical),
    unit_peak = max(abs(canonical)),
    unit_integral = sum(diff(ref_grid) *
      (utils::head(canonical, -1L) + utils::tail(canonical, -1L)) / 2),
    unit_peak_per_basis = if (is.matrix(ref_vals)) {
      apply(abs(ref_vals), 2L, max)
    } else {
      max(abs(canonical))
    }
  )

  if (identical(mode, "unit_peak_per_basis")) {
    factor[factor == 0] <- 1
    if (any(!is.finite(factor))) {
      stop("Per-basis normalization factors must be finite.", call. = FALSE)
    }
  } else if (!is.finite(factor) || abs(factor) < .Machine$double.xmin) {
    stop(
      sprintf("Normalization factor for '%s' (mode = '%s') is not usable.",
              orig_name, mode),
      call. = FALSE
    )
  }

  normalized_func <- function(t) {
    values <- hrf(t)
    if (length(factor) == 1L) {
      values / factor
    } else if (is.matrix(values)) {
      sweep(values, 2L, factor, "/")
    } else {
      values / factor
    }
  }

  as_hrf(
    f = normalized_func,
    name = paste0(orig_name, "[norm=", mode, "]"),
    nbasis = orig_nbasis,
    span = orig_span,
    params = c(decorator_params, list(
      .normalization = mode,
      .normalization_factor = factor
    ))
  )
}


#' Normalise Each Basis of an HRF to Unit Peak
#'
#' Back-compatible spelling and behavior for independently peak-normalizing
#' every basis column. New code can use
#' `normalize_hrf(hrf, "unit_peak_per_basis")` explicitly.
#'
#' @param hrf An object of class `HRF`.
#' @return A unit-peak `HRF` object.
#' @family HRF_decorator_functions
#' @export
#' @examples
#' gauss_unnorm <- as_hrf(function(t) 5 * dnorm(t, 6, 2), name = "unnorm_gauss")
#' gauss_norm <- normalise_hrf(gauss_unnorm)
#' max(gauss_norm(seq(0, 20, by = 0.1)))
normalise_hrf <- function(hrf) {
  out <- normalize_hrf(hrf, "unit_peak_per_basis")
  attr(out, "name") <- paste0(attr(hrf, "name"), "_norm")
  params <- attr(out, "params")
  params$.normalised <- TRUE
  attr(out, "params") <- params
  out
}
