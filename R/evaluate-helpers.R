#' @importFrom memoise memoise
#' @keywords internal
#' @noRd
# Bounded in-memory LRU cache for HRF evaluations to prevent
# unbounded growth in long-running sessions. If cachem isn't
# available at runtime, fall back to memoise's default memory cache.
.hrf_eval_cache <- tryCatch(
  cachem::cache_mem(
    max_size = getOption("fmrihrf.hrf_cache_max_size", 50 * 1024^2), # 50 MB default
    evict = "lru"
  ),
  error = function(e) memoise::cache_memory()
)

.memo_hrf <- memoise::memoise(
  function(hrf, span, dt) {
    if (!is.numeric(span) || length(span) != 1 || span <= 0) {
        stop("`span` must be a single numeric value strictly greater than 0.", call. = FALSE)
    }
    if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) {
        stop("`dt` must be a single numeric value strictly greater than 0.", call. = FALSE)
    }
    times <- seq(0, span, by = dt)
    # Evaluate HRF - ensure it returns a matrix
    val <- fmrihrf::evaluate(hrf, times)
    if (is.vector(val)) matrix(val, ncol = 1) else val
  },
  cache = .hrf_eval_cache
)

#' Prepare Inputs for Regressor Evaluation Engines
#'
#' Internal helper function to perform common setup steps before calling
#' a specific evaluation engine (fft, conv, loop, Rconv).
#' Handles filtering of events, evaluation/memoization of HRF on fine grid.
#'
#' @param x A `Reg` object.
#' @param grid The target evaluation time grid (numeric vector).
#' @param precision The precision for internal calculations (numeric scalar).
#' @return A list containing prepared inputs:
#'   * `nb`: Number of basis functions.
#'   * `hrf_span`: The span of the HRF.
#'   * `valid_ons`: Filtered onset times relevant to the grid.
#'   * `valid_durs`: Corresponding durations.
#'   * `valid_amp`: Corresponding amplitudes.
#'   * `grid`: The original target grid.
#'   * `precision`: The precision value.
#'   * `hrf_fine_matrix`: HRF values evaluated on the fine time grid (potentially memoized).
#'     NULL if hrf_is_list is TRUE.
#'   * `fine_grid`: The fine time grid itself (if needed by Rconv/loop).
#'   * `summate`: Logical summation flag from the regressor.
#'   * `hrf`: The original HRF object (single HRF or list of HRFs).
#'   * `hrf_is_list`: Logical indicating if hrf is a list of per-event HRFs.
#'   * `valid_hrfs`: If hrf_is_list, the subset of HRFs corresponding to valid events.
#' @keywords internal
#' @noRd
#' @importFrom stats approx median convolve
prep_reg_inputs <- function(x, grid, precision) {

  # Ensure grid is sorted (Correctness 1.4)
  if (is.unsorted(grid)) {
      warning("Input grid is unsorted. Sorting grid for evaluation.")
      grid <- sort(grid)
  }

  # Check if HRF is a list (trial-varying case)
  hrf_is_list <- isTRUE(attr(x, "hrf_is_list"))

  # Get nbasis - for list HRFs, use the first one (all should have same nbasis)
  if (hrf_is_list) {
    nb <- if (length(x$hrf) > 0) nbasis(x$hrf[[1]]) else 1L
  } else {
    nb <- nbasis(x$hrf)
  }

  hrf_span <- x$span

  # Filter events based on grid boundaries and HRF span
  onset_min_bound <- grid[1] - hrf_span
  onset_max_bound <- grid[length(grid)]

  # Start with potentially already filtered data from Reg constructor
  keep_indices <- which(x$onsets >= onset_min_bound & x$onsets <= onset_max_bound)

  # Note: Amplitude filtering already done in Reg(), no need to repeat here
  valid_ons <- x$onsets[keep_indices]
  valid_durs <- x$duration[keep_indices]
  valid_amp <- x$amplitude[keep_indices]

  # Also filter HRF list if applicable
  valid_hrfs <- if (hrf_is_list) x$hrf[keep_indices] else NULL

  if (length(valid_ons) == 0) {
    # Return minimal info needed to signal zero output
    return(list(nb = nb, grid = grid, valid_ons = numeric(0), hrf_is_list = hrf_is_list))
  }

  # Prepare/Memoize finely sampled HRF for efficient evaluation
  # Only for single HRF case - list HRFs are handled per-event in eval_loop
  if (hrf_is_list) {
    hrf_fine_matrix <- NULL
  } else {
    hrf_fine_matrix <- .memo_hrf(x$hrf, hrf_span, precision)
  }

  # Prepare fine grid (needed for Rconv/loop interpolation)
  # Start at the earliest possible contribution from kept events
  # which occurs at `grid[1] - hrf_span`
  fine_grid_start <- grid[1] - hrf_span
  # Use the full range of onsets when determining the end of the fine grid
  # to handle unsorted event inputs without reordering
  fine_grid_end <- max(grid[length(grid)], max(valid_ons) + max(valid_durs)) + hrf_span
  fine_grid <- seq(fine_grid_start, fine_grid_end, by = precision)

  return(list(
    nb         = nb,
    hrf_span   = hrf_span,
    valid_ons  = valid_ons,
    valid_durs = valid_durs,
    valid_amp  = valid_amp,
    grid       = grid,
    precision  = precision,
    hrf_fine_matrix = hrf_fine_matrix,
    fine_grid  = fine_grid,
    summate    = x$summate,
    hrf        = x$hrf,
    hrf_is_list = hrf_is_list,
    valid_hrfs = valid_hrfs
  ))
}

# Internal Evaluation Engines -----

#' FFT-based Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
eval_fft <- function(p, ...) {
  # List HRFs require per-event evaluation - fallback to loop

  if (isTRUE(p$hrf_is_list)) {
    return(eval_loop(p, ...))
  }

  # Call the unified C++ wrapper
  result <- evaluate_regressor_cpp(
              grid = p$grid,
              onsets = p$valid_ons,
              durations = p$valid_durs,
              amplitudes = p$valid_amp,
              hrf_matrix = p$hrf_fine_matrix,
              hrf_span = p$hrf_span,
              precision = p$precision,
              method = "fft"
            )
  result
}

#' Direct Convolution Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
eval_conv <- function(p, ...) {
  # List HRFs require per-event evaluation - fallback to loop
  if (isTRUE(p$hrf_is_list)) {
    return(eval_loop(p, ...))
  }

  # Call the unified C++ wrapper
  result <- evaluate_regressor_cpp(
              grid = p$grid,
              onsets = p$valid_ons,
              durations = p$valid_durs,
              amplitudes = p$valid_amp,
              hrf_matrix = p$hrf_fine_matrix,
              hrf_span = p$hrf_span,
              precision = p$precision,
              method = "conv"
            )
  result
}

#' R Convolution Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
#' @importFrom stats convolve approx
eval_Rconv <- function(p, ...) {
  # List HRFs require per-event evaluation - fallback to loop
  if (isTRUE(p$hrf_is_list)) {
    return(eval_loop(p, ...))
  }

  # Check conditions (moved from evaluate.Reg)
  is_regular_grid <- length(p$grid) > 1 && length(unique(round(diff(p$grid), 8))) == 1
  is_constant_duration <- length(unique(p$valid_durs)) <= 1

  if (!is_regular_grid || !is_constant_duration) {
    warning("Method 'Rconv' requires a regular grid and constant event durations. Falling back to 'loop' method.")
    return(eval_loop(p, ...)) # Call the loop engine directly as fallback
  }
  
  # Proceed with R convolution using stats::convolve
  delta <- numeric(length(p$fine_grid))
  onset_indices <- floor((p$valid_ons - p$fine_grid[1]) / p$precision) + 1
  valid_onset_indices <- onset_indices >= 1 & onset_indices <= length(p$fine_grid)

  if (length(p$valid_durs) > 0) {
    dur_len <- floor(p$valid_durs[1] / p$precision)
  } else {
    dur_len <- 0
  }

  for (i in which(valid_onset_indices)) {
    start_idx <- onset_indices[i]
    end_idx <- min(start_idx + dur_len, length(p$fine_grid))
    delta[start_idx:end_idx] <- delta[start_idx:end_idx] + p$valid_amp[i]
  }
  
  samhrf <- p$hrf_fine_matrix # Already evaluated and potentially memoized
  nb <- p$nb
  
  if (nb > 1) {
    lowres <- matrix(0, length(p$grid), nb)
    for (b in 1:nb) {
      highres_conv <- stats::convolve(delta, rev(samhrf[, b]), type = "open")
      valid_len <- length(p$fine_grid)
      highres_trimmed <- highres_conv[1:valid_len]
      interp_res <- approx(p$fine_grid, highres_trimmed, xout = p$grid, rule = 2)$y
      lowres[, b] <- interp_res
    }
    result <- lowres
  } else {
    highres_conv <- stats::convolve(delta, rev(as.vector(samhrf)), type = "open")
    valid_len <- length(p$fine_grid)
    highres_trimmed <- highres_conv[1:valid_len]
    result <- approx(p$fine_grid, highres_trimmed, xout = p$grid, rule = 2)$y
  }
  result
}

#' R Loop Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments passed to evaluate.HRF.
#' @keywords internal
#' @noRd
eval_loop <- function(p, ...) {
  # Check if HRF is a list (trial-varying case)
  hrf_is_list <- isTRUE(p$hrf_is_list)

  # Validate HRF(s)
  if (hrf_is_list) {
    if (is.null(p$valid_hrfs) || length(p$valid_hrfs) == 0) {
      stop("Error inside eval_loop: p$valid_hrfs is NULL or empty for list HRF case.")
    }
  } else {
    if (is.null(p$hrf) || !inherits(p$hrf, "HRF")) {
      stop("Error inside eval_loop: p$hrf is NULL or not an HRF object.")
    }
  }

  nb <- p$nb
  hrf_span <- p$hrf_span
  grid <- p$grid
  valid_ons <- p$valid_ons
  valid_durs <- p$valid_durs
  valid_amp <- p$valid_amp
  precision <- p$precision
  summate <- p$summate

  dspan <- hrf_span / stats::median(diff(grid), na.rm=TRUE) # Approx span in grid units

  # Pre-calculate nearest grid indices for onsets (more robust than RANN for this)
  # Find the index of the grid point *just before or at* each onset
  nidx <- findInterval(valid_ons, grid)
  nidx[nidx == 0] <- 1

  outmat <- matrix(0, length(grid), nb)

  for (i in seq_along(valid_ons)) {
    start_grid_idx <- nidx[i]
    end_grid_idx <- min(start_grid_idx + ceiling(dspan) + 5, length(grid))
    if (start_grid_idx > length(grid)) next
    grid.idx <- start_grid_idx:end_grid_idx

    relOns <- grid[grid.idx] - valid_ons[i]
    valid_rel_idx <- which(relOns >= 0 & relOns <= hrf_span)

    if (length(valid_rel_idx) > 0) {
        target_indices_outmat <- grid.idx[valid_rel_idx]

        # Select the appropriate HRF for this event
        current_hrf <- if (hrf_is_list) p$valid_hrfs[[i]] else p$hrf

        # Call evaluate S3 generic, should dispatch to evaluate.HRF
        resp <- evaluate(current_hrf, relOns[valid_rel_idx], amplitude=valid_amp[i],
                         duration=valid_durs[i],
                         precision=precision,
                         summate=summate, ...)

        if (!is.matrix(resp) && nb > 1) {
            resp <- matrix(resp, ncol=nb)
        }
        if (!is.matrix(resp) && nb == 1) {
            resp <- matrix(resp, ncol=1)
        }

        if (nrow(resp) != length(target_indices_outmat)){
            warning("Dimension mismatch between response and target indices in loop.")
            next
        }

        outmat[target_indices_outmat, seq_len(nb)] <-
          outmat[target_indices_outmat, seq_len(nb), drop = FALSE] + resp
    }
  }

  outmat
}
