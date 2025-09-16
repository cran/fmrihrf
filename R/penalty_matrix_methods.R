#' @rdname penalty_matrix
#' @export
penalty_matrix.HRF <- function(x, order = 2, ...) {
  diag(nbasis(x))
}

# internal helper for roughness penalties
roughness_penalty <- function(nb, order = 2) {
  if (nb <= 1) {
    diag(nb)
  } else if (nb > order) {
    D <- diff(diag(nb), differences = order)
    crossprod(D)
  } else {
    diag(nb)
  }
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.BSpline_HRF <- function(x, order = 2, ...) {
  roughness_penalty(nbasis(x), order)
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.Tent_HRF <- function(x, order = 2, ...) {
  roughness_penalty(nbasis(x), order)
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.FIR_HRF <- function(x, order = 2, ...) {
  roughness_penalty(nbasis(x), order)
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.SPMG2_HRF <- function(x, order = 2, shrink_deriv = 2, ...) {
  nb <- nbasis(x)
  R <- diag(nb)
  if (nb >= 1) R[1, 1] <- 0
  if (nb >= 2) R[2, 2] <- shrink_deriv
  R
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.SPMG3_HRF <- function(x, order = 2, shrink_deriv = 2, ...) {
  nb <- nbasis(x)
  R <- diag(nb)
  if (nb >= 1) R[1, 1] <- 0
  if (nb >= 2) R[2, 2] <- shrink_deriv
  if (nb >= 3) R[3, 3] <- shrink_deriv
  R
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.Fourier_HRF <- function(x, order = 2, ...) {
  nb <- nbasis(x)
  freqs <- ceiling(seq_len(nb) / 2)
  diag(freqs^order)
}

#' @rdname penalty_matrix
#' @export
penalty_matrix.Daguerre_HRF <- function(x, order = 2, ...) {
  nb <- nbasis(x)
  diag((seq_len(nb) - 1)^2)
}
