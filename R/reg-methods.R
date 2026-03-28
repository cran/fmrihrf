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
#' @param normalize Logical; if TRUE, scale evaluated regressor output to unit peak
#'   (maximum absolute value of 1). For multi-basis regressors, each basis column
#'   is normalized independently.
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
evaluate.Reg <- function(x, grid, precision=.33, method=c("conv", "fft", "Rconv", "loop"),
                         sparse = FALSE, normalize = FALSE, ...) {

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
  if (!is.logical(normalize) || length(normalize) != 1 || is.na(normalize)) {
    stop("`normalize` must be a single logical value.", call. = FALSE)
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

  if (normalize) {
    if (is.matrix(final_result)) {
      peaks <- apply(final_result, 2, function(col) max(abs(col), na.rm = TRUE))
      peaks[is.na(peaks) | peaks == 0] <- 1
      final_result <- sweep(final_result, 2, peaks, "/")
    } else {
      peak <- max(abs(final_result), na.rm = TRUE)
      if (!is.na(peak) && peak != 0) {
        final_result <- final_result / peak
      }
    }
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
#' @examples
#' r <- regressor(onsets = c(1, 10, 20), hrf = HRF_SPMG1,
#'                duration = 0, amplitude = 1,
#'                span = 40)
#' print(r)
#' @export
#' @method print Reg
#' @rdname print
print.Reg <- function(x, ...) {

  n_ons <- length(x$onsets)
  hrf_is_list <- isTRUE(attr(x, "hrf_is_list"))

  # Get HRF info - handle list vs single HRF

  if (hrf_is_list) {
    if (length(x$hrf) > 0) {
      hrf_name <- paste0("trial-varying (", length(x$hrf), " HRFs)")
      nb <- nbasis(x$hrf[[1]])
    } else {
      hrf_name <- "trial-varying (empty)"
      nb <- 1L
    }
  } else {
    hrf_name <- attr(x$hrf, "name") %||% "custom function"
    nb <- nbasis(x$hrf)
  }
  hrf_span <- x$span

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
nbasis.Reg <- function(x, ...) {
  # Handle list HRFs (trial-varying case)
  if (isTRUE(attr(x, "hrf_is_list"))) {
    if (length(x$hrf) > 0) {
      return(nbasis(x$hrf[[1]]))
    } else {
      return(1L)
    }
  }
  nbasis(x$hrf)
}

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

#' Plot a Regressor Object
#'
#' Creates a visualization of a regressor object showing the predicted BOLD
#' response over time. Optionally displays event onsets as vertical lines.
#'
#' @param x A `Reg` object created by `regressor()`.
#' @param grid Numeric vector of time points for evaluation. If NULL (default),
#'   automatically generates a grid from 0 to max(onsets) + span with step 0.5s.
#' @param show_onsets Logical; if TRUE (default), show vertical dashed lines
#'   at event onset times.
#' @param onset_color Color for onset lines. Default is "red".
#' @param onset_alpha Alpha transparency for onset lines. Default is 0.5.
#' @param precision Numeric sampling precision for HRF evaluation. Default is 0.33.
#' @param ... Additional arguments passed to underlying plot functions.
#' @return Invisibly returns a data frame with the time and response values.
#' @examples
#' # Create and plot a simple regressor
#' reg <- regressor(onsets = c(10, 30, 50), hrf = HRF_SPMG1)
#' plot(reg)
#'
#' # Plot with custom time grid
#' plot(reg, grid = seq(0, 80, by = 1))
#'
#' # Plot without onset markers
#' plot(reg, show_onsets = FALSE
#' )
#' @method plot Reg
#' @export
plot.Reg <- function(x, grid = NULL, show_onsets = TRUE,
                     onset_color = "red", onset_alpha = 0.5,
                     precision = 0.33, ...) {

  # Generate default grid if not provided
  if (is.null(grid)) {
    max_time <- max(x$onsets, na.rm = TRUE) + x$span
    grid <- seq(0, max_time, by = 0.5)
  }


  # Evaluate the regressor
  response <- evaluate(x, grid, precision = precision)


  # Get regressor info for title
  hrf_is_list <- isTRUE(attr(x, "hrf_is_list"))
  if (hrf_is_list) {
    hrf_name <- paste0("trial-varying (", length(x$hrf), " HRFs)")
  } else {
    hrf_name <- attr(x$hrf, "name") %||% "custom"
  }
  n_events <- length(x$onsets)
  title <- sprintf("Regressor: %d events, HRF: %s", n_events, hrf_name)

  if (is.matrix(response)) {
    # Multi-basis regressor
    nb <- ncol(response)
    graphics::matplot(grid, response, type = "l", lwd = 1.5, lty = 1,
                      xlab = "Time (s)", ylab = "Response",
                      main = title, ...)
    if (show_onsets) {
      graphics::abline(v = x$onsets, lty = 2, col = onset_color,
                       lwd = 0.5)
    }
    graphics::legend("topright", paste("Basis", 1:nb),
                     col = 1:nb, lty = 1, lwd = 1.5, bty = "n")

    df <- data.frame(time = grid, response)
    colnames(df)[-1] <- paste0("basis_", 1:nb)
  } else {
    # Single-basis regressor
    graphics::plot(grid, response, type = "l", lwd = 1.5,
                   xlab = "Time (s)", ylab = "Response",
                   main = title, ...)
    if (show_onsets) {
      graphics::abline(v = x$onsets, lty = 2, col = onset_color,
                       lwd = 0.5)
    }

    df <- data.frame(time = grid, response = response)
  }

  invisible(df)
}

#' Compare Multiple Regressor Objects
#'
#' Creates a comparison plot of multiple regressor objects. This function provides
#' a convenient way to visualize different regressors on the same plot, with options
#' for showing event onsets and customization. Uses ggplot2 if available for
#' publication-quality figures, otherwise falls back to base R graphics.
#'
#' @param ... Regressor objects to compare. Can be passed as individual arguments
#'   or as a named list.
#' @param grid Numeric vector of time points for evaluation. If NULL (default),
#'   automatically generates a grid covering all regressors.
#' @param labels Character vector of labels for each regressor. If NULL (default),
#'   uses "Regressor_1", "Regressor_2", etc.
#' @param title Character string for the plot title. If NULL (default),
#'   uses "Regressor Comparison".
#' @param subtitle Character string for the plot subtitle. If NULL, no subtitle.
#' @param show_onsets Logical or character. If TRUE, show onset lines for all
#'   regressors. If "first", show only for the first regressor. If FALSE, hide
#'   onsets. Default is "first".
#' @param onset_alpha Alpha transparency for onset lines. Default is 0.3.
#' @param precision Numeric sampling precision for HRF evaluation. Default is 0.33.
#' @param use_ggplot Logical; if TRUE and ggplot2 is available, use ggplot2
#'   for plotting. If FALSE, use base R graphics. Default is TRUE.
#' @return Invisibly returns a data frame in long format with columns 'time',
#'   'Regressor', and 'response'.
#' @examples
#' # Create regressors with different HRFs
#' onsets <- c(10, 30, 50)
#' reg1 <- regressor(onsets, HRF_SPMG1)
#' reg2 <- regressor(onsets, HRF_GAMMA)
#' reg3 <- regressor(onsets, HRF_GAUSSIAN)
#'
#' # Compare regressors
#' plot_regressors(reg1, reg2, reg3,
#'                 labels = c("SPM Canonical", "Gamma", "Gaussian"))
#'
#' # Compare regressors with different event timings
#' reg_fast <- regressor(seq(0, 60, by = 10), HRF_SPMG1)
#' reg_slow <- regressor(seq(0, 60, by = 20), HRF_SPMG1)
#' plot_regressors(reg_fast, reg_slow,
#'                 labels = c("Fast (10s ISI)", "Slow (20s ISI)"),
#'                 title = "Effect of Inter-Stimulus Interval")
#'
#' # Compare original vs shifted regressor
#' reg_orig <- regressor(c(10, 30, 50), HRF_SPMG1)
#' reg_shifted <- shift(reg_orig, 5)
#' plot_regressors(reg_orig, reg_shifted,
#'                 labels = c("Original", "Shifted +5s"))
#' @export
plot_regressors <- function(..., grid = NULL, labels = NULL,
                            title = NULL, subtitle = NULL,
                            show_onsets = "first", onset_alpha = 0.3,
                            precision = 0.33, use_ggplot = TRUE) {

  # Collect regressor objects
  regs <- list(...)

  # Handle case where a single list is passed
  if (length(regs) == 1 && is.list(regs[[1]]) && !inherits(regs[[1]], "Reg")) {
    regs <- regs[[1]]
  }

  n_regs <- length(regs)
  if (n_regs == 0) {
    stop("At least one Reg object must be provided")
  }

  # Validate all inputs are Reg objects
  for (i in seq_along(regs)) {
    if (!inherits(regs[[i]], "Reg")) {
      stop("All arguments must be Reg objects. Argument ", i, " is not a Reg.")
    }
  }

  # Determine time grid
  if (is.null(grid)) {
    # Find max time across all regressors
    max_times <- sapply(regs, function(r) {
      if (length(r$onsets) == 0) return(r$span)
      max(r$onsets, na.rm = TRUE) + r$span
    })
    max_time <- max(max_times)
    grid <- seq(0, max_time, by = 0.5)
  }

  # Generate labels
  if (is.null(labels)) {
    labels <- paste0("Regressor_", seq_along(regs))
  }
  if (length(labels) != n_regs) {
    stop("Length of 'labels' must match number of regressors")
  }

  # Evaluate all regressors
  responses <- lapply(regs, function(r) {
    resp <- evaluate(r, grid, precision = precision)
    # For multi-basis, take the first basis
    if (is.matrix(resp)) {
      resp[, 1]
    } else {
      resp
    }
  })

  # Build data frame
  df_list <- lapply(seq_along(responses), function(i) {
    data.frame(time = grid, Regressor = labels[i], response = responses[[i]])
  })
  df <- do.call(rbind, df_list)
  df$Regressor <- factor(df$Regressor, levels = labels)  # Preserve order

  # Set default title
  if (is.null(title)) {
    title <- "Regressor Comparison"
  }

  # Determine which onsets to show
  if (isTRUE(show_onsets)) {
    onset_data <- do.call(rbind, lapply(seq_along(regs), function(i) {
      data.frame(onset = regs[[i]]$onsets, Regressor = labels[i])
    }))
  } else if (identical(show_onsets, "first")) {
    onset_data <- data.frame(onset = regs[[1]]$onsets, Regressor = labels[1])
  } else {
    onset_data <- NULL
  }

  # Plot
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  if (use_ggplot && has_ggplot) {
    # ggplot2 version
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = response, color = Regressor)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = "Time (seconds)",
        y = "Response",
        color = "Regressor"
      ) +
      ggplot2::theme_minimal()

    # Add onset lines if requested
    if (!is.null(onset_data) && nrow(onset_data) > 0) {
      p <- p + ggplot2::geom_vline(
        data = onset_data,
        ggplot2::aes(xintercept = onset),
        linetype = "dashed", alpha = onset_alpha, color = "gray40"
      )
    }

    print(p)
    attr(df, "plot") <- p
  } else {
    # Base R version
    colors <- grDevices::rainbow(n_regs)
    y_range <- range(df$response, na.rm = TRUE)

    graphics::plot(NULL, xlim = range(grid), ylim = y_range,
                   xlab = "Time (seconds)", ylab = "Response",
                   main = title)
    if (!is.null(subtitle)) {
      graphics::mtext(subtitle, side = 3, line = 0.5, cex = 0.8)
    }

    # Add onset lines first (so they're behind the data)
    if (!is.null(onset_data) && nrow(onset_data) > 0) {
      graphics::abline(v = unique(onset_data$onset), lty = 2,
                       col = grDevices::adjustcolor("gray40", alpha.f = onset_alpha))
    }

    for (i in seq_along(responses)) {
      graphics::lines(grid, responses[[i]], col = colors[i], lwd = 1.5)
    }

    graphics::legend("topright", legend = labels, col = colors, lty = 1,
                     lwd = 1.5, bty = "n")
  }

  invisible(df)
}
