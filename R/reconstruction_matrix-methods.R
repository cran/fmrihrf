#' Reconstruction matrix for an HRF object
#'
#' @description
#' S3 method for `HRF` objects that returns a matrix mapping basis coefficients
#' to sampled HRF values at the provided time grid. For single-basis HRFs, this
#' returns a one-column matrix. For multi-basis HRFs (e.g., SPMG2/SPMG3, FIR,
#' B-spline), this returns a matrix with one column per basis function.
#'
#' @param hrf An object of class `HRF`.
#' @param sframe A numeric vector of times, or a `sampling_frame` object from which
#'   times are extracted via `samples()`.
#' @param ... Additional arguments passed to `samples()` when `sframe` is a
#'   `sampling_frame`, and to `evaluate()` for HRF evaluation.
#'
#' @return A numeric matrix of dimension `length(times) x nbasis(hrf)`.
#' @rdname reconstruction_matrix
#' @method reconstruction_matrix HRF
#' @export
reconstruction_matrix.HRF <- function(hrf, sframe, ...) {
  # Derive time grid from numeric vector or sampling_frame
  times <- if (inherits(sframe, "sampling_frame")) {
    samples(sframe, ...)
  } else if (is.numeric(sframe)) {
    sframe
  } else {
    stop("`sframe` must be a numeric vector of times or a `sampling_frame` object.", call. = FALSE)
  }

  vals <- evaluate(hrf, times, ...)
  if (is.matrix(vals)) {
    vals
  } else {
    matrix(vals, ncol = 1)
  }
}

