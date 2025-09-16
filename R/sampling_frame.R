#' @noRd
#' @keywords internal
new_sampling_frame <- function(blocklens, TR, start_time, precision) {
  structure(
    list(blocklens = blocklens,
         TR        = TR,
         start_time = start_time,
         precision = precision),
    class = "sampling_frame")
}

#'
#' A \code{sampling_frame} describes the block structure and temporal sampling of an fMRI paradigm.
#'
#' @param blocklens A numeric vector representing the number of scans in each block.
#' @param TR A numeric value or vector representing the repetition time in seconds
#'   (i.e., the spacing between consecutive image acquisitions). When a vector is
#'   provided, its length must be 1 or equal to the number of blocks.
#' @param start_time A numeric value or vector representing the offset of the first
#'   scan of each block (default is \code{TR/2}). When a vector is provided, its
#'   length must be 1 or equal to the number of blocks.
#' @param precision A numeric value representing the discrete sampling interval used for convolution with the hemodynamic response function (default is 0.1).
#'
#' @examples
#' frame <- sampling_frame(blocklens = c(100, 100, 100), TR = 2, precision = 0.5)
#'
#' # The relative time (with respect to the last block) in seconds of each sample/acquisition
#' sam <- samples(frame)
#' # The global time (with respect to the first block) of each sample/acquisition
#' gsam <- samples(frame, global = TRUE)
#'
#' # Block identifiers for each acquisition can be retrieved using
#' # blockids(frame)
#'
#' @return A list with class "sampling_frame" describing the block structure and temporal sampling of an fMRI paradigm.
#' @export
sampling_frame <- function(blocklens, TR, start_time = TR / 2, precision = .1)
{
  # --- recycle & validate ------------------------------------------------
  n_blocks <- length(blocklens)

  # Check that TR and start_time are either length 1 or match n_blocks
  if (length(TR) != 1L && length(TR) != n_blocks) {
    stop("TR must have length 1 or match the number of blocks")
  }
  if (length(start_time) != 1L && length(start_time) != n_blocks) {
    stop("start_time must have length 1 or match the number of blocks")
  }

  # Recycle values to the number of blocks
  TR <- rep_len(TR, n_blocks)
  start_time <- rep_len(start_time, n_blocks)
  
  # Validate inputs with proper error messages
  if (!is.numeric(blocklens) || any(is.na(blocklens))) {
    stop("Block lengths must be numeric and non-NA")
  }
  if (!is.numeric(TR) || any(is.na(TR))) {
    stop("TR values must be numeric and non-NA")
  }
  if (!is.numeric(start_time) || any(is.na(start_time))) {
    stop("Start times must be numeric and non-NA")
  }
  if (!all(blocklens > 0)) {
    stop("Block lengths must be positive")
  }
  if (!all(TR > 0)) {
    stop("TR values must be positive")
  }
  if (!all(start_time >= 0)) {
    stop("Start times must be non-negative")
  }
  if (precision <= 0) {
    stop("Precision must be positive")
  }
  if (precision >= min(TR)) {
    stop("Precision must be positive and less than the minimum TR")
  }

  new_sampling_frame(blocklens, TR, start_time, precision)
}


#' @method samples sampling_frame
#' @rdname samples
#' @export
samples.sampling_frame <- function(x, blockids = NULL, global = FALSE,...) {
  if (is.null(blockids)) blockids <- seq_along(x$blocklens)

  # number of scans per selected block
  lens <- x$blocklens[blockids]

  # fast allocate
  idx <- sequence(lens) - 1                   # 0‑based within block
  
  # Calculate relative times within each block
  block_times <- rep(blockids, lens)  # which block each sample belongs to
  times <- x$start_time[block_times] + idx * x$TR[block_times]

  if (global) {
    # For global timing, add the cumulative time offset from previous blocks
    # Calculate cumulative time at the end of each block
    block_durations <- x$blocklens * x$TR
    cumulative_time <- c(0, cumsum(block_durations))
    
    # Add the cumulative time offset for each block
    time_offsets <- cumulative_time[block_times]
    times + time_offsets
  } else {
    times
  }
}

#' @method global_onsets sampling_frame
#' @rdname global_onsets
#' @param blockids Integer vector identifying the block for each onset. Values
#'   must be whole numbers with no NAs.
#' @export
global_onsets.sampling_frame <- function(x, onsets, blockids,...) {
  # Calculate cumulative time offsets for each block
  block_durations <- x$blocklens * x$TR
  cumulative_time <- c(0, cumsum(block_durations))

  if (!is.numeric(blockids) || anyNA(blockids) || any(blockids %% 1 != 0)) {
    stop("blockids must be whole numbers and not NA")
  }
  blockids <- as.integer(blockids)
  stopifnot(length(onsets) == length(blockids),
            all(blockids >= 1L), all(blockids <= length(x$blocklens)))

  onsets + cumulative_time[blockids]
}

#' @export
#' @rdname print
print.sampling_frame <- function(x, ...) {
  n_blk <- length(x$blocklens)
  total_scans <- sum(x$blocklens)
  
  cat("Sampling Frame\n")
  cat("==============\n\n")
  
  cat("Structure:\n")
  cat(sprintf("  %d block%s\n", n_blk, if (n_blk > 1) "s" else ""))
  cat(sprintf("  Total scans: %d\n\n", total_scans))
  
  cat("Timing:\n")
  cat(sprintf("  TR: %s s\n", paste(unique(x$TR), collapse = ", ")))
  cat(sprintf("  Precision: %.3g s\n\n", x$precision))
  
  cat("Duration:\n")
  total_time <- sum(x$blocklens * x$TR)
  cat(sprintf("  Total time: %.1f s\n", total_time))
  
  invisible(x)
}


#' @rdname blockids
#' @export
blockids.sampling_frame <- function(x, ...) {
  rep(seq_along(x$blocklens), times = x$blocklens)
}


#' @rdname blocklens
#' @export
blocklens.sampling_frame <- function(x,...) {
    x$blocklens
}


#' @rdname acquisition_onsets
#' @method acquisition_onsets sampling_frame
#' @export
acquisition_onsets.sampling_frame <- function(x, ...) {
  samples(x, global = TRUE)
}
