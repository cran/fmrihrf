#' Default derivative method for HRF objects
#'
#' Uses numerical differentiation via numDeriv::grad when analytic derivatives
#' are not available for a specific HRF type.
#'
#' @param x An HRF object
#' @param t Numeric vector of time points at which to evaluate the derivative
#' @param ... Additional arguments (currently unused)
#' @return Numeric vector or matrix of derivative values
#' @importFrom numDeriv grad
#' @examples
#' t <- seq(0, 30, by = 0.5)
#' d <- deriv(HRF_SPMG1, t)
#' @export
#' @method deriv HRF
deriv.HRF <- function(x, t, ...) {
  # Get the number of basis functions
  nb <- nbasis(x)
  
  # Create result matrix
  result <- matrix(0, nrow = length(t), ncol = nb)
  
  # For multi-basis HRFs, compute derivative for each basis function separately
  if (nb > 1) {
    for (j in 1:nb) {
      # Define function for j-th basis function
      f_j <- function(time) {
        val <- evaluate(x, time)
        if (is.matrix(val)) {
          return(val[1, j])
        } else {
          return(val[1])
        }
      }
      
      # Compute derivative for each time point
      for (i in seq_along(t)) {
        result[i, j] <- numDeriv::grad(f_j, t[i])
      }
    }
    return(result)
  } else {
    # For single basis functions, use simpler approach
    for (i in seq_along(t)) {
      t_i <- t[i]
      
      # Define function to differentiate at time t_i
      f <- function(time) {
        val <- evaluate(x, time)
        if (is.matrix(val)) {
          return(val[1, 1])
        } else {
          return(val[1])
        }
      }
      
      # Compute gradient at t_i
      result[i, 1] <- numDeriv::grad(f, t_i)
    }
    
    return(result[, 1])
  }
}


#' Derivative method for SPMG1 HRF
#'
#' Uses the analytic derivative formula for the SPM canonical HRF.
#'
#' @param x An SPMG1_HRF object
#' @param t Numeric vector of time points at which to evaluate the derivative
#' @param ... Additional arguments (currently unused)
#' @return Numeric vector of derivative values
#' @rdname deriv.HRF
#' @export
#' @method deriv SPMG1_HRF
deriv.SPMG1_HRF <- function(x, t, ...) {
  # Extract parameters from the HRF object
  params <- attr(x, "params")
  if (is.null(params)) {
    # Use defaults if params not found
    params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  }
  
  # Call the analytic derivative function
  hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
}


#' Derivative method for SPMG2 HRF
#'
#' Returns derivatives for both the canonical HRF and its temporal derivative.
#' The first column contains the derivative of the canonical HRF, and the second
#' column contains the second derivative (derivative of the temporal derivative).
#'
#' @param x An SPMG2_HRF object
#' @param t Numeric vector of time points at which to evaluate the derivative
#' @param ... Additional arguments (currently unused)
#' @return Matrix with 2 columns of derivative values
#' @rdname deriv.HRF
#' @export
#' @method deriv SPMG2_HRF
deriv.SPMG2_HRF <- function(x, t, ...) {
  # Extract parameters - SPMG2 is a basis set, so get params from first basis
  params <- NULL
  basis_list <- attr(x, "basis_list")
  if (!is.null(basis_list) && length(basis_list) >= 1) {
    params <- attr(basis_list[[1]], "params")
  }
  
  if (is.null(params)) {
    # Use defaults if params not found
    params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  }
  
  # Create result matrix with 2 columns
  result <- matrix(0, nrow = length(t), ncol = 2)
  
  # First column: derivative of canonical HRF
  result[, 1] <- hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  
  # Second column: second derivative (derivative of temporal derivative)
  result[, 2] <- hrf_spmg1_second_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  
  return(result)
}


#' Derivative method for SPMG3 HRF
#'
#' Returns derivatives for the canonical HRF and its two derivatives.
#' Since SPMG3 already includes first and second derivatives as basis functions,
#' this method returns their derivatives (1st, 2nd, and 3rd derivatives of the original HRF).
#'
#' @param x An SPMG3_HRF object
#' @param t Numeric vector of time points at which to evaluate the derivative
#' @param ... Additional arguments (currently unused)
#' @return Matrix with 3 columns of derivative values
#' @rdname deriv.HRF
#' @export
#' @method deriv SPMG3_HRF
deriv.SPMG3_HRF <- function(x, t, ...) {
  # Extract parameters - SPMG3 is a basis set, so get params from first basis
  params <- NULL
  basis_list <- attr(x, "basis_list")
  if (!is.null(basis_list) && length(basis_list) >= 1) {
    params <- attr(basis_list[[1]], "params")
  }
  
  if (is.null(params)) {
    # Use defaults if params not found
    params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  }
  
  # Create result matrix with 3 columns
  result <- matrix(0, nrow = length(t), ncol = 3)
  
  # First column: derivative of canonical HRF (1st derivative)
  result[, 1] <- hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  
  # Second column: second derivative
  result[, 2] <- hrf_spmg1_second_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  
  # Third column: third derivative (use numerical differentiation of second derivative)
  # Define function for second derivative
  f_second_deriv <- function(time) {
    hrf_spmg1_second_deriv(time, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  }
  
  # Compute third derivative numerically
  for (i in seq_along(t)) {
    result[i, 3] <- numDeriv::grad(f_second_deriv, t[i])
  }
  
  return(result)
}