# Evaluate method for Reg objects
#' 
#' This is the primary method for evaluating regressor objects created by the `Reg` constructor
#' (and thus also works for objects created by `regressor`).
#' It dispatches to different internal methods based on the `method` argument.
#' 
#' @rdname evaluate
#' @param x A `Reg` object (or an object inheriting from it, like `regressor`).
#' @param grid Numeric vector specifying the time points (seconds) for evaluation.
#' @param precision Numeric sampling precision for internal HRF evaluation and convolution (seconds).
#' @param method The evaluation method: 
#'   \describe{
#'     \item{conv}{(Default) Uses the C++ direct convolution (`evaluate_regressor_convolution`). Generally safer and more predictable.}
#'     \item{fft}{Uses the fast C++ FFT convolution (`evaluate_regressor_fast`). Can be faster but may fail with very fine precision or wide grids.  
#'       Extremely fine `precision` or wide `grid` ranges may trigger an internal FFT size exceeding ~1e7, which results in an error.}
#'     \item{Rconv}{Uses an R-based convolution (`stats::convolve`). Requires constant event durations and a regular sampling grid. Can be faster than the R loop for many events meeting these criteria.}
#'     \item{loop}{Uses a pure R implementation involving looping through onsets. Can be slower, especially for many onsets.}
#'   }
#' @param sparse Logical indicating whether to return a sparse matrix (from the Matrix package). Default is FALSE.
#' @param ... Additional arguments passed down (e.g., to `evaluate.HRF` in the loop method).
#' @examples
#' # Create a regressor
#' reg <- regressor(onsets = c(10, 30, 50), hrf = HRF_SPMG1)
#' 
#' # Evaluate with default method (conv)
#' times <- seq(0, 80, by = 0.5)
#' response <- evaluate(reg, times)
#' 
#' # Try different evaluation methods
#' response_loop <- evaluate(reg, times, method = "loop")
#' 
#' # With higher precision
#' response_precise <- evaluate(reg, times, precision = 0.1)
#' @export
#' @method evaluate Reg
#' @importFrom Matrix Matrix
#' @importFrom memoise memoise
#' @importFrom stats approx median convolve
#' @importFrom Rcpp evalCpp
evaluate.Reg <- function(x, grid, precision=.33, method=c("conv", "fft", "Rconv", "loop"), sparse = FALSE, ...) {

  method <- match.arg(method)

  # Validate inputs
  if (!is.numeric(grid) || length(grid) == 0 || anyNA(grid)) {
    stop("`grid` must be a non-empty numeric vector with no NA values.",
         call. = FALSE)
  }
  if (!is.numeric(precision) || length(precision) != 1 || is.na(precision) ||
      precision <= 0) {
    stop("`precision` must be a positive numeric value.", call. = FALSE)
  }

  # Prepare inputs using the helper function
  prep_data <- prep_reg_inputs(x, grid, precision)
  
  # Check if prep_reg_inputs indicated no relevant events
  if (length(prep_data$valid_ons) == 0) {
    zero_res <- if (prep_data$nb == 1) {
      rep(0, length(grid))
    } else {
      matrix(0, nrow = length(grid), ncol = prep_data$nb)
    }
    if (sparse) {
      zero_res <- Matrix::Matrix(zero_res, sparse = TRUE)
    }
    return(zero_res)
  }
  
  # --- Method Dispatch to Internal Engines ---
  eng_fun <- switch(method,
     conv  = eval_conv,   # Now the default - safer direct convolution
     fft   = eval_fft,    # FFT-based (faster but can fail with large FFT sizes)
     Rconv = eval_Rconv,  # R-based convolution
     loop  = eval_loop,   # Pure R loop implementation
     stop("Invalid evaluation method: ", method) # Should not happen due to match.arg
  )
  
  # Call the selected engine function with prepared data
  # Pass ... through to the engine, which might pass it to evaluate.HRF in loop
  result <- eng_fun(prep_data, ...) 
  
  # --- Final Formatting ---
  nb <- prep_data$nb
  final_result <- if (nb == 1 && is.matrix(result)) {
    as.vector(result)
  } else if (nb > 1 && !is.matrix(result)) {
    matrix(result, nrow=length(grid), ncol=nb)
  } else {
      result
  }
  
  # Convert to sparse matrix if requested
  if (sparse) {
    return(Matrix::Matrix(final_result, sparse = TRUE))
  } else {
    return(final_result)
  }
}


#' @method shift Reg
#' @rdname shift
#' @export
#' @importFrom assertthat assert_that
shift.Reg <- function(x, shift_amount, ...) {
  dots <- list(...)

  if (missing(shift_amount) && "offset" %in% names(dots)) {
    shift_amount <- dots$offset
  }

  assert_that(inherits(x, "Reg"),
              msg = "Input 'x' must inherit from class 'Reg'")

  if (missing(shift_amount)) {
    stop("Must supply `shift_amount` or `offset`.", call. = FALSE)
  }

  assert_that(is.numeric(shift_amount) && length(shift_amount) == 1,
              msg = "`shift_amount` must be a single numeric value")

  # Handle empty regressor case
  if (length(x$onsets) == 0 || (length(x$onsets) == 1 && is.na(x$onsets[1]))) {
    # Returning the original empty object is appropriate for a shift
    return(x)
  }

  # Shift the valid onsets
  shifted_onsets <- x$onsets + shift_amount

  # Reconstruct the object using the core Reg constructor 
  out <- Reg(onsets = shifted_onsets,
             hrf = x$hrf,
             duration = x$duration,
             amplitude = x$amplitude,
             span = x$span,
             summate = x$summate)
             
  return(out)
}

#' Print method for Reg objects
#' 
#' Provides a concise summary of the regressor object using the cli package.
#' 
#' @param x A `Reg` object.
#' @param ... Not used.
#' @return No return value, called for side effects (prints to console)
#' @importFrom cli cli_h1 cli_text cli_div cli_li
#' @export
#' @method print Reg
#' @rdname print
print.Reg <- function(x, ...) {
  
  n_ons <- length(x$onsets)
  hrf_name <- attr(x$hrf, "name") %||% "custom function"
  nb <- nbasis(x$hrf)
  hrf_span <- attr(x$hrf, "span") %||% x$span
  
  cli::cli_h1("fMRI Regressor Object")
  
  # Use cli_div for potentially better alignment than cli_ul
  cli::cli_div(theme = list(ul = list("margin-left" = 2), li = list("margin-bottom" = 0.5)))
  cli::cli_li("Type: {.cls {class(x)[1]}}{if(inherits(x, 'regressor')) ' (Legacy compatible)'}")
  if (n_ons == 0) {
    cli::cli_li("Events: 0 (Empty Regressor)")
  } else {
    cli::cli_li("Events: {n_ons}")
    cli::cli_li("Onset Range: {round(min(x$onsets), 2)}s to {round(max(x$onsets), 2)}s")
    if (any(x$duration != 0)) {
      cli::cli_li("Duration Range: {round(min(x$duration), 2)}s to {round(max(x$duration), 2)}s")
    }
    if (!all(x$amplitude == 1)) {
      cli::cli_li("Amplitude Range: {round(min(x$amplitude), 2)} to {round(max(x$amplitude), 2)}")
    }
  }
  cli::cli_li("HRF: {hrf_name} ({nb} basis function{?s})")
  cli::cli_li("HRF Span: {hrf_span}s")
  cli::cli_li("Summation: {x$summate}")
  cli::cli_end()

  invisible(x)
}

# S3 Methods for Reg class -----

#' @export
#' @rdname nbasis
#' @method nbasis Reg
nbasis.Reg <- function(x, ...) nbasis(x$hrf)

#' @export
#' @rdname onsets
#' @method onsets Reg
onsets.Reg <- function(x, ...) x$onsets

#' @export
#' @rdname durations
#' @method durations Reg
durations.Reg <- function(x, ...) x$duration

#' @export
#' @rdname amplitudes
#' @method amplitudes Reg
amplitudes.Reg <- function(x, ...) x$amplitude
