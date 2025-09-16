#' Construct a Regressor Set
#'
#' Creates a set of regressors, one for each level of a factor. Each
#' condition shares the same HRF and other parameters but has distinct
#' onsets, durations and amplitudes.
#'
#' @param onsets Numeric vector of event onset times.
#' @param fac A factor (or object coercible to a factor) indicating the
#'   condition for each onset.
#' @param hrf Hemodynamic response function used for all conditions.
#' @param duration Numeric scalar or vector of event durations.
#' @param amplitude Numeric scalar or vector of event amplitudes.
#' @param span Numeric scalar giving the HRF span in seconds.
#' @param summate Logical; passed to [regressor()].
#'
#' @return An object of class `RegSet` containing one `Reg` per factor level.
#' @examples
#' # Create events for 3 conditions
#' onsets <- c(10, 20, 30, 40, 50, 60)
#' conditions <- factor(c("A", "B", "C", "A", "B", "C"))
#' 
#' # Create regressor set
#' rset <- regressor_set(onsets, conditions, hrf = HRF_SPMG1)
#' 
#' # With durations and amplitudes
#' rset2 <- regressor_set(
#'   onsets = onsets,
#'   fac = conditions,
#'   duration = 2,
#'   amplitude = c(1, 1.5, 0.8, 1, 1.5, 0.8),
#'   hrf = HRF_SPMG1
#' )
#' 
#' # Evaluate the regressor set
#' times <- seq(0, 80, by = 0.1)
#' design_matrix <- evaluate(rset, times)
#' @export
regressor_set <- function(onsets, fac, hrf = HRF_SPMG1, duration = 0,
                          amplitude = 1, span = 40, summate = TRUE) {
  fac <- as.factor(fac)
  onsets    <- as.numeric(onsets)
  duration  <- recycle_or_error(as.numeric(duration), length(onsets), "duration")
  amplitude <- recycle_or_error(as.numeric(amplitude), length(onsets), "amplitude")
  if (length(fac) != length(onsets)) {
    stop("`fac` must be the same length as `onsets`", call. = FALSE)
  }
  levs <- levels(fac)
  regs <- lapply(levs, function(lv) {
    idx <- which(fac == lv)
    if (length(idx) == 0) {
      null_regressor(hrf = hrf, span = span)
    } else {
      Reg(onsets = onsets[idx], hrf = hrf,
          duration = duration[idx], amplitude = amplitude[idx],
          span = span, summate = summate)
    }
  })
  structure(list(regs = regs, levels = levs), class = c("RegSet", "list"))
}

#' @rdname regressor_set
#' @param x A RegSet object
#' @param grid Numeric vector of time points at which to evaluate
#' @param precision Numeric precision for evaluation
#' @param method Evaluation method
#' @param sparse Logical whether to return sparse matrix
#' @param ... Additional arguments passed to evaluate
#' @export
evaluate.RegSet <- function(x, grid, precision = .33,
                            method = c("conv", "fft", "Rconv", "loop"),
                            sparse = FALSE, ...) {
  method <- match.arg(method)
  mats <- lapply(x$regs, evaluate, grid = grid, precision = precision,
                 method = method, sparse = FALSE, ...)
  out <- do.call(cbind, mats)
  if (sparse) {
    return(Matrix::Matrix(out, sparse = TRUE))
  } else {
    return(out)
  }
}
