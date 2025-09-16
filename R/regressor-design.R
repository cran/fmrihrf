#' Build a Design Matrix from Block-wise Onsets
#'
#' `regressor_design` extends [regressor_set()] by allowing onsets to be
#' specified relative to individual blocks and by directly returning the
#' evaluated design matrix.
#'
#' @param onsets Numeric vector of event onset times, expressed relative to the
#'   start of their corresponding block.
#' @param fac A factor (or object coercible to a factor) indicating the
#'   condition for each onset.
#' @param block Integer vector identifying the block for each onset. Values must
#'   be valid block indices for `sframe`.
#' @param sframe A [sampling_frame] describing the temporal structure of the
#'   experiment.
#' @param hrf Hemodynamic response function shared by all conditions.
#' @param duration Numeric scalar or vector of event durations.
#' @param amplitude Numeric scalar or vector of event amplitudes.
#' @param span Numeric scalar giving the HRF span in seconds.
#' @param precision Numeric precision used during convolution.
#' @param method Evaluation method passed to [evaluate()].
#' @param sparse Logical; if `TRUE` a sparse design matrix is returned.
#' @param summate Logical; passed to [regressor()].
#'
#' @return A numeric matrix (or sparse matrix) with one column per factor level
#'   and one row per sample defined by `sframe`.
#' @examples
#' # Create a sampling frame for 2 blocks, 100 scans each, TR=2
#' sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
#' 
#' # Events in block-relative time
#' onsets <- c(10, 30, 50, 20, 40, 60)
#' conditions <- factor(c("A", "B", "A", "B", "A", "B"))
#' blocks <- c(1, 1, 1, 2, 2, 2)
#' 
#' # Build design matrix
#' design <- regressor_design(
#'   onsets = onsets,
#'   fac = conditions,
#'   block = blocks,
#'   sframe = sframe,
#'   hrf = HRF_SPMG1
#' )
#' 
#' # Design matrix has 200 rows (total scans) and 2 columns (conditions)
#' dim(design)
#' @export
regressor_design <- function(onsets, fac, block, sframe,
                             hrf = HRF_SPMG1, duration = 0,
                             amplitude = 1, span = 40,
                             precision = .33,
                             method = c("conv", "fft", "Rconv", "loop"),
                             sparse = FALSE,
                             summate = TRUE) {
  fac   <- as.factor(fac)
  block <- as.integer(block)

  if (length(onsets) != length(fac) || length(onsets) != length(block)) {
    stop("`onsets`, `fac` and `block` must have the same length", call. = FALSE)
  }

  method <- match.arg(method)

  # Convert block-wise onsets to global timing
  g_onsets <- global_onsets(sframe, onsets, block)

  # Create regressor set using global onsets
  rs <- regressor_set(g_onsets, fac, hrf = hrf,
                      duration = duration, amplitude = amplitude,
                      span = span, summate = summate)

  # Evaluation grid corresponds to all acquisition times
  grid <- samples(sframe, global = TRUE)

  evaluate(rs, grid = grid, precision = precision,
           method = method, sparse = sparse)
}

