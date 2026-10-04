# Put this in R/utils-internal.R

# Declare global variables to avoid R CMD check NOTEs for ggplot2 aes() usage
utils::globalVariables(c("time", "response", "HRF", "Regressor", "onset", ".data"))

#' @keywords keyword
#' @noRd
recycle_or_error <- function(x, n, name) {
  if (length(x) == n) {
    return(x)
  }
  if (length(x) == 1) {
    return(rep(x, n))
  }
  stop(paste0("`", name, "` must have length 1 or ", n, ", not ", length(x)), call. = FALSE)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# Zero a response outside the causal support t >= 0.
#
# Several shape functions (gaussian, mexican hat, inverse logit, LWU) are
# defined by formulas that are non-zero at negative lag. The regressor path
# masks negative lags itself, but decorators such as block_hrf() and lag_hrf()
# sample h(t - offset) directly and would otherwise mix in pre-onset values.
#' @keywords internal
#' @noRd
.causal <- function(t, values) {
  values[!is.na(t) & t < 0] <- 0
  values
}

#' @importFrom utils tail
#' @keywords internal
#' @noRd
.block_offsets_weights <- function(width, precision) {
  offsets <- seq(0, width, by = precision)
  if (tail(offsets, 1) < width) {
    offsets <- c(offsets, width)
  }
  if (length(offsets) == 1) {
    return(list(offsets = offsets, weights = 1))
  }

  deltas <- diff(offsets)
  weights <- numeric(length(offsets))
  weights[1] <- deltas[1] / 2
  weights[length(offsets)] <- deltas[length(deltas)] / 2
  if (length(offsets) > 2) {
    weights[2:(length(offsets) - 1)] <- (deltas[-length(deltas)] + deltas[-1]) / 2
  }
  list(offsets = offsets, weights = weights)
}

# Combine a list of HRF evaluations (one per quadrature offset) into a single
# weighted result. Handles both single-basis (vector) and multi-basis (matrix)
# outputs uniformly. When `summate = FALSE`, the result is divided by the sum
# of weights so amplitude does not scale with block width.
#
# @param value_list list of numeric vectors (single-basis) or matrices
#   (multi-basis). All elements must share the same shape.
# @param weights numeric vector of quadrature weights, same length as value_list.
# @param nbasis optional declared number of basis functions. When supplied and
#   equal to 1, the result is flattened to a numeric vector (matching the
#   original single-basis dispatch in block_hrf). When NULL (default), dispatch
#   is based solely on the shape of value_list[[1]].
# @param summate logical; if FALSE, divide by sum(weights).
# @return a numeric vector or matrix matching the shape of value_list[[1]].
#' @keywords internal
#' @noRd
.weighted_combine <- function(value_list, weights, nbasis = NULL, summate = TRUE) {
  as_vector <- if (!is.null(nbasis)) {
    isTRUE(nbasis == 1L)
  } else {
    !is.matrix(value_list[[1]])
  }
  if (as_vector) {
    hmat <- do.call(cbind, lapply(value_list, as.numeric))
    res <- as.vector(hmat %*% weights)
  } else {
    weighted <- Map(function(vals, wt) vals * wt, value_list, weights)
    res <- Reduce("+", weighted)
  }
  if (!summate) {
    weight_sum <- sum(weights)
    if (weight_sum > 0) {
      res <- res / weight_sum
    }
  }
  res
}

# Peak-normalise a vector or matrix. For matrices, each column is normalised
# independently. Columns (or vectors) whose peak absolute value is 0 or NA are
# returned unchanged.
#' @keywords internal
#' @noRd
.normalise_result <- function(x) {
  if (is.matrix(x)) {
    peaks <- apply(x, 2, function(col) max(abs(col), na.rm = TRUE))
    peaks[is.na(peaks) | peaks == 0] <- 1
    sweep(x, 2, peaks, "/")
  } else {
    peak_val <- max(abs(x), na.rm = TRUE)
    if (!is.na(peak_val) && peak_val != 0) x / peak_val else x
  }
}

# Compute per-basis peak absolute values from a reference evaluation. Zero or
# NA peaks are replaced with 1 so downstream division is a no-op. Returns a
# vector of length ncol(x) for matrices, or a scalar for vectors.
#' @keywords internal
#' @noRd
.get_peaks <- function(x) {
  if (is.matrix(x)) {
    peaks <- apply(x, 2, function(col) max(abs(col), na.rm = TRUE))
  } else {
    peaks <- max(abs(as.numeric(x)), na.rm = TRUE)
  }
  peaks[is.na(peaks) | peaks == 0] <- 1
  peaks
}
