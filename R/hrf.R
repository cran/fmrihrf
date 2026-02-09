#' Turn any function into an HRF object
#'
#' This is the core constructor for creating HRF objects in the refactored system.
#' It takes a function `f(t)` and attaches standard HRF attributes. If `params` are
#' provided, `as_hrf` creates a closure that captures these parameters, ensuring they
#' are used during evaluation rather than relying on the function's defaults.
#'
#' @param f The function to be turned into an HRF object. It must accept a single argument `t` (time).
#' @param name The name for the HRF object. Defaults to the deparsed name of `f`.
#' @param nbasis The number of basis functions represented by `f`. Must be \code{>= 1}. Defaults to 1L.
#' @param span The nominal time span (duration in seconds) of the HRF. Must be positive. Defaults to 24.
#' @param params A named list of parameters associated with the HRF function `f`. When provided,
#'   `as_hrf` creates a closure that captures these parameters. Defaults to an empty list.
#' @return A new HRF object.
#' @examples
#' # Create a custom HRF from a function
#' custom_hrf <- as_hrf(function(t) exp(-t/5),
#'                      name = "exponential",
#'                      span = 20)
#' evaluate(custom_hrf, seq(0, 10, by = 1))
#'
#' # Create HRF with specific parameters (closure captures them)
#' gamma_hrf <- as_hrf(hrf_gamma, params = list(shape = 8, rate = 1.2))
#' evaluate(gamma_hrf, seq(0, 20, by = 1))
#' @keywords internal
#' @export
as_hrf <- function(f, name = deparse(substitute(f)), nbasis = 1L, span = 24,
                   params = list()) {
  assertthat::assert_that(is.function(f))
  assertthat::assert_that(is.character(name), length(name) == 1)
  assertthat::assert_that(is.numeric(nbasis), length(nbasis) == 1)
  assertthat::assert_that(nbasis >= 1, msg = "nbasis must be >= 1")
  assertthat::assert_that(is.numeric(span), length(span) == 1)
  assertthat::assert_that(span > 0, msg = "span must be > 0")
  assertthat::assert_that(is.list(params))

  # If params are provided, create a closure that captures them
  if (length(params) > 0) {
    # Validate that f accepts the provided parameters
    f_formals <- names(formals(f))
    param_names <- names(params)

    # Params starting with '.' are metadata (e.g., .lag, .width from decorators)
    # and should not be validated as function arguments
    metadata_params <- param_names[startsWith(param_names, ".")]
    callable_params <- param_names[!startsWith(param_names, ".")]

    # Check if callable param names are valid arguments to f (excluding 't')
    # 't' is special and should not be in params
    invalid_params <- setdiff(callable_params, f_formals)
    if (length(invalid_params) > 0) {
      warning(sprintf("Parameters %s are not arguments to function %s and will be ignored",
                      paste(invalid_params, collapse = ", "), name),
              call. = FALSE)
      # Keep metadata params, filter out invalid callable params
      params <- params[param_names %in% c(f_formals, metadata_params)]
    }

    # Separate callable params (to be passed to function) from metadata (just stored)
    callable_params_list <- params[!startsWith(names(params), ".")]

    # Create closure that captures parameters if any callable params remain
    if (length(callable_params_list) > 0) {
      # Capture the original function in a local variable to avoid recursion issues
      orig_f <- f
      # Create a closure that calls the original function with captured callable parameters only
      f <- function(t) {
        do.call(orig_f, c(list(t = t), callable_params_list))
      }
    }
  }

  structure(
    f,
    class        = c("HRF", "function"),
    name         = name,
    nbasis       = as.integer(nbasis),
    span         = span,
    param_names  = names(params),
    params       = params
  )
}


#' Bind HRFs into a Basis Set
#'
#' Combines multiple HRF objects into a single multi-basis HRF object.
#' The resulting function evaluates each input HRF at time `t` and returns the results column-bound together.
#'
#' @param ... One or more HRF objects created by `as_hrf` or other HRF constructors/decorators.
#'
#' @return A new HRF object representing the combined basis set.
#'
#' @examples
#' # Combine multiple HRF basis functions
#' hrf1 <- as_hrf(hrf_gaussian, params = list(mean = 5))
#' hrf2 <- as_hrf(hrf_gaussian, params = list(mean = 10))
#' basis <- bind_basis(hrf1, hrf2)
#' nbasis(basis)  # Returns 2
#' 
#' @keywords internal
#' @export
#' @importFrom assertthat assert_that
bind_basis <- function(...) {
  xs <- list(...)
  assertthat::assert_that(length(xs) > 0, msg = "bind_basis requires at least one HRF object.")
  assertthat::assert_that(all(sapply(xs, inherits, "HRF")), msg = "All inputs to bind_basis must be HRF objects.")

  # Handle single HRF case explicitly
  if (length(xs) == 1) {
    return(xs[[1]])
  }

  # Calculate combined attributes
  combined_nbasis <- sum(vapply(xs, attr, 0L, "nbasis"))
  combined_span <- max(vapply(xs, attr, 0, "span"))
  combined_name <- paste(sapply(xs, attr, "name"), collapse = " + ")

  # Create the combined function
  combined_func <- function(t) {
    do.call(cbind, lapply(xs, function(f) f(t)))
  }

  # Use as_hrf to create the new HRF object
  as_hrf(
    f = combined_func,
    name = combined_name,
    nbasis = combined_nbasis,
    span = combined_span,
    params = list() # Params usually don't combine meaningfully
  )
}


#' Construct an HRF Instance using Decorators
#' 
#' @description
#' `gen_hrf` takes a base HRF function or object and applies optional lag,
#' blocking, and normalization decorators based on arguments.
#'
#' @param hrf A function `f(t)` or an existing `HRF` object.
#' @param lag Optional lag in seconds. If non-zero, applies `lag_hrf`.
#' @param width Optional block width in seconds. If non-zero, applies `block_hrf`.
#' @param precision Sampling precision for block convolution (passed to `block_hrf`). Default is 0.1.
#' @param half_life Half-life decay parameter for exponential decay in seconds (passed to `block_hrf`). Default is Inf (no decay).
#' @param summate Whether to summate within blocks (passed to `block_hrf`). Default is TRUE.
#' @param normalize If TRUE, applies `normalise_hrf` at the end. Default is FALSE.
#' @param name Optional name for the *final* HRF object. If NULL (default), a name is generated based on the base HRF and applied decorators.
#' @param span Optional span for the *final* HRF object. If NULL (default), the span is determined by the base HRF and decorators.
#' @param ... Extra arguments passed to the *base* HRF function if `hrf` is a function.
#'
#' @return A final `HRF` object, potentially modified by decorators.
#' 
#' @examples 
#' # Lagged SPMG1
#' grf_lag <- gen_hrf(HRF_SPMG1, lag=3)
#' # Blocked Gaussian
#' grf_block <- gen_hrf(hrf_gaussian, width=5, precision=0.2)
#' # Lagged and Blocked, then Normalized
#' grf_both_norm <- gen_hrf(HRF_SPMG1, lag=2, width=4, normalize=TRUE)
#'
#' @export
gen_hrf <- function(hrf, lag=0, width=0, precision=.1, half_life=Inf,
                    summate=TRUE, normalize=FALSE, name=NULL, span=NULL, ...) {

  # 1. Ensure we start with an HRF object
  if (is.function(hrf) && !inherits(hrf, "HRF")) {
    # If it's a plain function, convert it using as_hrf
    # Determine nbasis by evaluating the function
    test_t <- 1:10 # A small sample range
    test_val <- try(hrf(test_t, ...), silent = TRUE)
    determined_nbasis <- if (!inherits(test_val, "try-error") && !is.null(test_val)) {
      if (is.matrix(test_val)) ncol(test_val) else 1L
    } else {
      warning(paste("Could not determine nbasis for function", deparse(substitute(hrf)), "- defaulting to 1. Evaluation failed."))
      1L
    }
    
    # Pass extra args (...) here if they are meant for the base function construction
    base_hrf <- as_hrf(f = function(t) hrf(t, ...),
                       name = deparse(substitute(hrf)),
                       nbasis = determined_nbasis) # Pass determined nbasis
                       # Let as_hrf determine default span, params
  } else if (inherits(hrf, "HRF")) {
    # If already an HRF object, use it directly
    base_hrf <- hrf
    if (length(list(...)) > 0) {
      warning("Ignoring extra arguments (...) because 'hrf' is already an HRF object.")
    }
  } else {
    stop("'hrf' must be a function or an HRF object.")
  }

  # Apply decorators conditionally
  decorated_hrf <- base_hrf

  # Apply width decorator first if needed
  if (width != 0) {
    # Check positivity *before* applying
    stopifnot(width > 0)
    # Note: block_hrf handles normalize=FALSE internally by default
    decorated_hrf <- block_hrf(decorated_hrf, width=width, precision=precision,
                               half_life=half_life, summate=summate, normalize=FALSE)
  }

  # Apply lag decorator if needed
  if (lag != 0) {
    decorated_hrf <- lag_hrf(decorated_hrf, lag=lag)
  }

  # Apply normalization decorator last if needed
  if (normalize) {
    decorated_hrf <- normalise_hrf(decorated_hrf)
  }

  # Override name and span if provided by user
  if (!is.null(name)) {
    attr(decorated_hrf, "name") <- name
  }
  if (!is.null(span)) {
    attr(decorated_hrf, "span") <- span
  }

  # Return the final (potentially decorated) HRF object
  return(decorated_hrf)
}


#' Generate an Empirical Hemodynamic Response Function
#' 
#' @description
#' `empirical_hrf` generates an empirical HRF using provided time points and values.
#' 
#' @param t Time points.
#' @param y Values of HRF at time `t[i]`.
#' @param name Name of the generated HRF.
#' @return An instance of type `HRF`.
#' @examples
#' # Create empirical HRF from data points
#' t_points <- seq(0, 20, by = 1)
#' y_values <- c(0, 0.1, 0.5, 0.9, 1.0, 0.8, 0.5, 0.2, 0, -0.1, -0.1, 
#'               0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' emp_hrf <- empirical_hrf(t_points, y_values)
#' 
#' # Evaluate at new time points
#' new_times <- seq(0, 25, by = 0.1)
#' response <- evaluate(emp_hrf, new_times)
#' @export
empirical_hrf <- function(t, y, name = "empirical_hrf") {
  as_hrf(stats::approxfun(t, y, yright = 0, yleft = 0),
         name = name, nbasis = 1L, span = max(t, na.rm = TRUE))
}

#' @export
#' @rdname empirical_hrf
#' @keywords internal
gen_empirical_hrf <- function(...) {
  .Deprecated("empirical_hrf")
  empirical_hrf(...)
}


#' Generate an HRF Basis Set
#' 
#' @description
#' `hrf_set` constructs an HRF basis set from one or more component HRF objects.
#'
#' @param ... One or more HRF objects.
#' @param name The name for the combined HRF set.
#' @return A combined HRF object.
#' @examples
#' # Combine multiple HRF types into a basis set
#' hrf_basis <- hrf_set(HRF_SPMG1, HRF_GAUSSIAN, HRF_GAMMA)
#' 
#' # Create custom basis with different parameters
#' hrf1 <- gen_hrf(hrf_gamma, alpha = 6, beta = 1)
#' hrf2 <- gen_hrf(hrf_gamma, alpha = 8, beta = 1)
#' custom_basis <- hrf_set(hrf1, hrf2, name = "custom_gamma_basis")
#' 
#' # Evaluate the basis set
#' t <- seq(0, 30, by = 0.1)
#' basis_response <- evaluate(hrf_basis, t)
#' @export
hrf_set <- function(..., name = "hrf_set") {
  combined_hrf <- bind_basis(...)
  attr(combined_hrf, "name") <- name
  combined_hrf
}

#' @export
#' @rdname hrf_set
#' @keywords internal
gen_hrf_set <- function(...) {
  .Deprecated("hrf_set")
  hrf_set(...)
}


#' Generate an HRF library from a parameter grid
#'
#' @description
#' `hrf_library` applies a base HRF generating function to each row of a parameter grid.
#'
#' @param fun A function that generates an HRF, given a set of parameters.
#' @param pgrid A data frame where each row is a set of parameters.
#' @param ... Additional arguments passed to `fun`.
#' @return A combined HRF object representing the library.
#' @examples
#' # Create library of gamma HRFs with varying parameters
#' param_grid <- expand.grid(
#'   shape = c(6, 8, 10),
#'   rate = c(0.9, 1, 1.1)
#' )
#' gamma_library <- hrf_library(
#'   function(shape, rate) as_hrf(hrf_gamma, params = list(shape = shape, rate = rate)),
#'   param_grid
#' )
#' 
#' # Create library with fixed and varying parameters
#' param_grid2 <- expand.grid(lag = c(0, 2, 4))
#' lagged_library <- hrf_library(
#'   function(lag) gen_hrf(HRF_SPMG1, lag = lag),
#'   param_grid2
#' )
#' @importFrom purrr pmap partial
#' @export
hrf_library <- function(fun, pgrid, ...) {
  extras <- list(...)
  # Ensure fun returns an HRF object
  hrf_list <- purrr::pmap(pgrid, function(...) {
      params <- list(...)
      do.call(fun, c(params, extras))
  })
  # Bind the generated HRFs
  do.call(bind_basis, hrf_list)
}

#' @export
#' @rdname hrf_library
#' @keywords internal
gen_hrf_library <- function(...) {
  .Deprecated("hrf_library")
  hrf_library(...)
}


#' HRF Constructor Function
#'
#' The `HRF` function creates an object representing a hemodynamic response function (HRF). It is a class constructor for HRFs.
#'
#' @param fun A function representing the hemodynamic response, mapping from time to BOLD response.
#' @param name A string specifying the name of the function.
#' @param nbasis An integer representing the number of basis functions, e.g., the columnar dimension of the HRF. Default is 1.
#' @param span A numeric value representing the span in seconds of the HRF. Default is 24.
#' @param param_names A character vector containing the names of the parameters for the HRF function.
#'
#' @return An HRF object with the specified properties.
#'
#' @details
#' The package provides several pre-defined HRF types that can be used in modeling fMRI responses:
#'
#' **Canonical HRFs:**
#' * `"spmg1"` or `HRF_SPMG1`: SPM's canonical HRF (single basis function)
#' * `"spmg2"` or `HRF_SPMG2`: SPM canonical + temporal derivative (2 basis functions)
#' * `"spmg3"` or `HRF_SPMG3`: SPM canonical + temporal and dispersion derivatives (3 basis functions)
#' * `"gaussian"` or `HRF_GAUSSIAN`: Gaussian-shaped HRF with peak around 5-6s
#' * `"gamma"` or `HRF_GAMMA`: Gamma function-based HRF with longer tail
#'
#' **Flexible basis sets:**
#' * `"bspline"` or `"bs"` or `HRF_BSPLINE`: B-spline basis for flexible HRF modeling
#' * `"tent"`: Tent (triangular) basis functions for flexible HRF modeling
#' * `"daguerre"` or `HRF_DAGUERRE`: Daguerre basis functions
#'
#' To see a complete list of available HRF types with details, use the `list_available_hrfs()` function.
#'
#' @examples
#' hrf <- HRF(hrf_gamma, "gamma", nbasis=1, param_names=c("shape", "rate"))
#' resp <- evaluate(hrf, seq(0, 24, by=1))
#'
#' # List all available HRF types
#' list_available_hrfs(details = TRUE)
#'
#' @export
#' @rdname HRF-class
HRF <- function(fun, name, nbasis=1, span=24, param_names=NULL) {
  vals <- try(fun(seq(0, span)), silent = TRUE)

  peak <- if (!inherits(vals, "try-error") && !is.null(vals)) {
    if (nbasis == 1) {
      max(vals, na.rm = TRUE)
    } else if (is.matrix(vals)) {
      max(apply(vals, 2, max, na.rm = TRUE))
    } else {
      NA # Unable to determine peak
    }
  } else {
    NA # Error during evaluation or null result
  }

  scale_factor <- if (!is.na(peak) && peak != 0) {
    1 / peak
  } else {
    NA # Cannot compute scale_factor if peak is NA or zero
  }
  
  structure(fun, name=name, 
            nbasis=as.integer(nbasis), 
            span=span,
            param_names=param_names, 
            scale_factor=scale_factor, 
            class=c("HRF", "function"))
  
}

#' @rdname nbasis
#' @export
nbasis.HRF <- function(x,...) attr(x, "nbasis")



#' @keywords internal
#' @noRd
makeDeriv <- function(HRF, n=1) {
  #with_package("numDeriv")
  if (n == 1) {
    function(t) numDeriv::grad(HRF, t)
  } else {
    Recall(function(t) numDeriv::grad(HRF,t), n-1)
  }
}


#' Generate a Lagged HRF Function
#'
#' @description
#' The `gen_hrf_lagged` function takes an HRF function and applies a specified lag to it. This can be useful for modeling time-delayed hemodynamic responses.
#'
#' @param hrf A function representing the underlying HRF to be shifted.
#' @param lag A numeric value specifying the lag or delay in seconds to apply to the HRF. This can also be a vector of lags, in which case the function returns an HRF set.
#' @param normalize A logical value indicating whether to rescale the output so that the maximum absolute value is 1. Defaults to `FALSE`.
#' @param ... Extra arguments supplied to the `hrf` function.
#'
#' @return A function representing the lagged HRF. If `lag` is a vector of lags, the function returns an HRF set.
#' @family gen_hrf
#' @examples
#' \donttest{
#' hrf_lag5 <- gen_hrf_lagged(HRF_SPMG1, lag=5)
#' hrf_lag5(0:20)
#' }
#'
#' @export
gen_hrf_lagged <- function(hrf, lag=2, normalize=FALSE, ...) {
  force(hrf)
  # TODO deal with nbasis arg in ...
  if (length(lag)>1) {
    do.call(gen_hrf_set, lapply(lag, function(l) gen_hrf_lagged(hrf, l,...)))
  } else {
    function(t) {
      ret <- hrf(t-lag,...)
      if (normalize) {
        ret <- ret/max(abs(ret))
      } 
      
      ret
    }
  }
}

#' @export
#' @describeIn gen_hrf_lagged alias for gen_hrf_lagged
#' @family gen_hrf
#' @return an lagged hrf function
hrf_lagged <- gen_hrf_lagged


#' Generate a Blocked HRF Function
#'
#' @description
#' The `gen_hrf_blocked` function creates a blocked HRF by convolving the input HRF with a boxcar function. This can be used to model block designs in fMRI analysis.
#'
#' @param hrf A function representing the hemodynamic response function. Default is `hrf_gaussian`.
#' @param width A numeric value specifying the width of the block in seconds. Default is 5.
#' @param precision A numeric value specifying the sampling resolution in seconds. Default is 0.1.
#' @param half_life A numeric value specifying the half-life of the exponential decay function, used to model response attenuation. Default is `Inf`, which means no decay.
#' @param summate A logical value indicating whether to allow each impulse response function to "add" up. Default is `TRUE`.
#' @param normalize A logical value indicating whether to rescale the output so that the peak of the output is 1. Default is `FALSE`.
#' @param ... Extra arguments passed to the HRF function.
#' @family gen_hrf
#'
#' @return A \code{function} representing the blocked HRF.
#' @examples
#' # Deprecated: use gen_hrf(..., width = 10) or block_hrf(HRF, width = 10)
#' @importFrom purrr partial
#' @export
gen_hrf_blocked <- function(hrf=hrf_gaussian, width=5, precision=.1, 
                            half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
  # Deprecated - use gen_hrf with width parameter instead
  .Deprecated("gen_hrf")
  gen_hrf(hrf, width=width, precision=precision, half_life=half_life, 
          summate=summate, normalize=normalize, ...)
}

#' @export
#' @aliases gen_hrf_blocked
#' @describeIn gen_hrf_blocked alias for gen_hrf_blocked
#' @return A \code{function} representing the blocked HRF.
hrf_blocked <- gen_hrf_blocked





#' Soft-threshold function
#'
#' This function applies soft-thresholding to the input values, setting values below the threshold to zero
#' and shrinking the remaining values by the threshold amount.
#'
#' @param x A numeric vector of input values
#' @param threshold A non-negative threshold value for the soft-thresholding operation
#'
#' @return A numeric vector with the soft-thresholded values
#'
#' @noRd
#' @keywords internal
soft_threshold <- function(x, threshold) {
  if (threshold < 0) {
    stop("Threshold value should be non-negative.")
  }

  sign(x) * pmax(0, abs(x) - threshold)
}



#' List all available hemodynamic response functions (HRFs)
#'
#' @description
#' Reads the internal HRF registry to list available HRF types.
#'
#' @param details Logical; if TRUE, attempt to add descriptions (basic for now).
#' @return A data frame with columns: name, type (object/generator), nbasis_default.
#' @examples
#' # List all available HRFs
#' hrfs <- list_available_hrfs()
#' print(hrfs)
#' 
#' # List with details
#' hrfs_detailed <- list_available_hrfs(details = TRUE)
#' print(hrfs_detailed)
#' @export
list_available_hrfs <- function(details = FALSE) {
  # Get names directly from the registry
  hrf_names <- names(HRF_REGISTRY)
  
  # Determine type and default nbasis by inspecting registry entries
  hrf_info <- lapply(hrf_names, function(name) {
    entry <- HRF_REGISTRY[[name]]
    type <- if (inherits(entry, 'HRF')) "object" else if (is.function(entry)) "generator" else "unknown"
    
    nbasis_default <- NA
    if (type == "object") {
      nbasis_default <- tryCatch(nbasis(entry), error = function(e) NA)
    } else if (type == "generator") {
      fmls <- formals(entry)
      if ("nbasis" %in% names(fmls)) {
        nb_val <- fmls$nbasis
        if(is.numeric(nb_val)) nbasis_default <- nb_val 
      } 
      if(is.na(nbasis_default)) nbasis_default <- "variable"
    }
    
    # Check if this name is an alias (points to the same object/func as another primary name)
    is_alias <- FALSE
    if (type == "object") {
      primary_names <- names(HRF_REGISTRY)[sapply(HRF_REGISTRY, identical, entry)]
      is_alias <- length(primary_names) > 1 && name %in% primary_names[primary_names != name]
    } else if (type == "generator") {
      # More complex for functions, check if it points to the same generator function
      # For now, let's assume aliases only exist for objects, or mark known ones
      is_alias <- name %in% c("gam", "bs")
    }

    list(name = name, type = type, nbasis_default = as.character(nbasis_default), is_alias = is_alias) 
  })
  
  # Combine into a data frame
  hrf_df <- do.call(rbind.data.frame, c(hrf_info, list(stringsAsFactors = FALSE)))
  
  # Add basic descriptions if requested
  if (details) {
      hrf_df$description <- paste(hrf_df$name, "HRF", 
                                  ifelse(hrf_df$type == "generator", "(generator)", "(object)"),
                                  ifelse(hrf_df$is_alias, "(alias)", ""))
  }
  
  hrf_df
}

# Define Static HRF Objects -----

#' Pre-defined Hemodynamic Response Function Objects
#' 
#' A collection of pre-defined HRF objects for common fMRI analysis scenarios.
#' These objects can be used directly in model specifications or as templates
#' for creating custom HRFs.
#' 
#' @section Canonical HRFs:
#' \describe{
#'   \item{\code{HRF_SPMG1}}{SPM canonical HRF (single basis function)}
#'   \item{\code{HRF_SPMG2}}{SPM canonical HRF with temporal derivative (2 basis functions)}
#'   \item{\code{HRF_SPMG3}}{SPM canonical HRF with temporal and dispersion derivatives (3 basis functions)}
#'   \item{\code{HRF_GAMMA}}{Gamma function-based HRF}
#'   \item{\code{HRF_GAUSSIAN}}{Gaussian function-based HRF}
#' }
#' 
#' @section Flexible Basis Sets:
#' \describe{
#'   \item{\code{HRF_BSPLINE}}{B-spline basis HRF (5 basis functions)}
#'   \item{\code{HRF_FIR}}{Finite Impulse Response (FIR) basis HRF (12 basis functions)}
#' }
#' 
#' @section Creating Custom Basis Sets:
#' The pre-defined objects above have fixed numbers of basis functions. To create
#' basis sets with custom parameters (e.g., different numbers of basis functions),
#' use one of these approaches:
#' 
#' \strong{Using getHRF():}
#' \itemize{
#'   \item \code{getHRF("fir", nbasis = 20)} - FIR basis with 20 functions
#'   \item \code{getHRF("bspline", nbasis = 10, span = 30)} - B-spline with 10 functions
#'   \item \code{getHRF("fourier", nbasis = 7)} - Fourier basis with 7 functions
#'   \item \code{getHRF("daguerre", nbasis = 5, scale = 3)} - Daguerre basis
#' }
#' 
#' \strong{Using generator functions directly:}
#' \itemize{
#'   \item \code{hrf_fir_generator(nbasis = 20, span = 30)}
#'   \item \code{hrf_bspline_generator(nbasis = 10, span = 30)}
#'   \item \code{hrf_fourier_generator(nbasis = 7, span = 24)}
#'   \item \code{hrf_daguerre_generator(nbasis = 5, scale = 3)}
#' }
#' 
#' @section Usage:
#' All HRF objects can be:
#' \itemize{
#'   \item Called as functions with time argument: \code{HRF_SPMG1(t)}
#'   \item Used in model specifications: \code{hrf(condition, basis = HRF_SPMG1)}
#'   \item Evaluated with \code{evaluate()} method
#'   \item Combined with decorators like \code{lag_hrf()} or \code{block_hrf()}
#' }
#' 
#' @param t Numeric vector of time points (in seconds) at which to evaluate the HRF
#'
#' @return 
#' When called as functions, return numeric vectors or matrices of HRF values.
#' When used as objects, they are HRF objects with class \code{c("HRF", "function")}.
#' 
#' @examples
#' # Evaluate HRFs at specific time points
#' times <- seq(0, 20, by = 0.5)
#' 
#' # Single basis canonical HRF
#' canonical_response <- HRF_SPMG1(times)
#' plot(times, canonical_response, type = "l", main = "SPM Canonical HRF")
#' 
#' # Multi-basis HRF with derivatives
#' multi_response <- HRF_SPMG3(times)  # Returns 3-column matrix
#' matplot(times, multi_response, type = "l", main = "SPM HRF with Derivatives")
#' 
#' # Gamma and Gaussian HRFs
#' gamma_response <- HRF_GAMMA(times)
#' gaussian_response <- HRF_GAUSSIAN(times)
#' 
#' # Compare different HRF shapes
#' plot(times, canonical_response, type = "l", col = "blue", 
#'      main = "HRF Comparison", ylab = "Response")
#' lines(times, gamma_response, col = "red")
#' lines(times, gaussian_response, col = "green")
#' legend("topright", c("SPM Canonical", "Gamma", "Gaussian"), 
#'        col = c("blue", "red", "green"), lty = 1)
#' 
#' # Create custom FIR basis with 20 bins
#' custom_fir <- getHRF("fir", nbasis = 20, span = 30)
#' fir_response <- evaluate(custom_fir, times)
#' matplot(times, fir_response, type = "l", main = "Custom FIR with 20 bins")
#' 
#' # Create custom B-spline basis  
#' custom_bspline <- hrf_bspline_generator(nbasis = 8, span = 25)
#' bspline_response <- evaluate(custom_bspline, times)
#' matplot(times, bspline_response, type = "l", main = "Custom B-spline with 8 basis functions")
#' 
#' @name HRF_objects
#' @aliases HRF_SPMG1 HRF_SPMG2 HRF_SPMG3 HRF_GAMMA HRF_GAUSSIAN HRF_BSPLINE HRF_FIR
#' @family hrf
#' @seealso 
#' \code{\link{evaluate.HRF}} for evaluating HRF objects,
#' \code{\link{gen_hrf}} for creating HRFs with decorators,
#' \code{\link{list_available_hrfs}} for listing all HRF types,
#' \code{\link{getHRF}} for creating HRFs by name with custom parameters,
#' \code{\link{hrf_fir_generator}}, \code{\link{hrf_bspline_generator}}, 
#' \code{\link{hrf_fourier_generator}}, \code{\link{hrf_daguerre_generator}} 
#' for creating custom basis sets directly
NULL

#' @rdname HRF_objects
#' @export
HRF_GAMMA <- as_hrf(hrf_gamma, name="gamma", params=list(shape=6, rate=1))

#' @rdname HRF_objects
#' @export
HRF_GAUSSIAN <- as_hrf(hrf_gaussian, name="gaussian", params=list(mean=6, sd=2))

#' @rdname HRF_objects
#' @export
HRF_SPMG1 <- as_hrf(hrf_spmg1, name="SPMG1", params=list(P1=5, P2=15, A1=0.0833))

#' @rdname HRF_objects
#' @export
HRF_SPMG2 <- bind_basis(
  as_hrf(hrf_spmg1, name="SPMG1_canonical", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_deriv, name="SPMG1_temporal_deriv", params=list(P1=5, P2=15, A1=0.0833))
)
attr(HRF_SPMG2, "name") <- "SPMG2"
class(HRF_SPMG2) <- c("SPMG2_HRF", class(HRF_SPMG2))

#' @rdname HRF_objects
#' @export
HRF_SPMG3 <- bind_basis(
  as_hrf(hrf_spmg1, name="SPMG1_canonical", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_deriv, name="SPMG1_temporal_deriv", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_second_deriv, name="SPMG1_dispersion_deriv", params=list(P1=5, P2=15, A1=0.0833))
)
attr(HRF_SPMG3, "name") <- "SPMG3"
class(HRF_SPMG3) <- c("SPMG3_HRF", class(HRF_SPMG3))

# Define HRF Generators (Functions returning HRF objects) -----

#' Create B-spline HRF Basis Set
#'
#' Generates an HRF object using B-spline basis functions with custom parameters.
#' This is the generator function that creates HRF objects with variable numbers
#' of basis functions, unlike the pre-defined \code{HRF_BSPLINE} which has 5 functions.
#'
#' @param nbasis Number of basis functions (default: 5)
#' @param span Temporal window in seconds (default: 24)
#' @return An HRF object of class \code{c("BSpline_HRF", "HRF", "function")}
#' @seealso \code{\link{HRF_objects}} for pre-defined HRF objects,
#'   \code{\link{getHRF}} for a unified interface to create HRFs
#' @examples
#' # Create B-spline basis with 10 functions
#' custom_bs <- hrf_bspline_generator(nbasis = 10)
#' t <- seq(0, 24, by = 0.1)
#' response <- evaluate(custom_bs, t)
#' matplot(t, response, type = "l", main = "B-spline HRF with 10 basis functions")
#' @export
hrf_bspline_generator <- function(nbasis=5, span=24) {
  # Validate inputs
  if (nbasis < 1) {
    stop("nbasis must be at least 1", call. = FALSE)
  }
  if (span <= 0) {
    stop("span must be positive", call. = FALSE)
  }
  
  # Ensure nbasis is integer
  nbasis <- as.integer(nbasis)
  
  degree <- 3 # Default cubic B-splines
  effective_nbasis <- max(1, nbasis) 
  
  f_bspline <- function(t) {
    valid_t_idx <- t >= 0 & t <= span
    if (!any(valid_t_idx)) {
      return(matrix(0, nrow = length(t), ncol = effective_nbasis))
    }
    
    res_mat <- matrix(0, nrow = length(t), ncol = effective_nbasis)
    
    bs_matrix <- tryCatch({
        splines::bs(t[valid_t_idx], df = effective_nbasis, degree = degree, 
                    Boundary.knots = c(0, span), intercept = FALSE)
    }, error = function(e) {
        warning(sprintf("splines::bs failed for effective_nbasis=%d, span=%d: %s", 
                        effective_nbasis, span, e$message), call. = FALSE)
        NULL 
    })

    if (!is.null(bs_matrix) && ncol(bs_matrix) == effective_nbasis) {
      res_mat[valid_t_idx, ] <- bs_matrix
    } else if (!is.null(bs_matrix)) {
       warning(sprintf("splines::bs returned %d columns, expected %d for effective_nbasis=%d, span=%d. Returning zeros.", 
                      ncol(bs_matrix), effective_nbasis, effective_nbasis, span), call. = FALSE)
    }
    return(res_mat)
  }

  obj <- as_hrf(
    f = f_bspline,
    name = "bspline", nbasis = as.integer(effective_nbasis), span = span
  )
  # Store params as metadata (the closure already captured these values)
  attr(obj, "params") <- list(nbasis = effective_nbasis, degree = degree, span = span)
  attr(obj, "param_names") <- c("nbasis", "degree", "span")
  class(obj) <- c("BSpline_HRF", class(obj))
  obj
}

#' Create Tent HRF Basis Set
#'
#' Generates an HRF object using tent (piecewise linear) basis functions with
#' custom parameters. This generator mirrors \code{HRF_TENT} but allows callers
#' to control the number of basis elements and temporal span.
#'
#' @param nbasis Number of tent basis functions (default: 5)
#' @param span Temporal window in seconds (default: 24)
#' @return An HRF object of class \code{c("Tent_HRF", "HRF", "function")}
#' @seealso \code{\link{HRF_objects}} for pre-defined HRF objects,
#'   \code{\link{getHRF}} for a unified interface to create HRFs,
#'   \code{\link{hrf_bspline_generator}} for a smoother alternative
#' @examples
#' # Create a tent basis with 6 functions over a 20 second window
#' custom_tent <- hrf_tent_generator(nbasis = 6, span = 20)
#' t <- seq(0, 20, by = 0.1)
#' response <- evaluate(custom_tent, t)
#' matplot(t, response, type = "l", main = "Tent HRF with 6 basis functions")
#' @export
hrf_tent_generator <- function(nbasis=5, span=24) {
  obj <- as_hrf(
    f = function(t) hrf_bspline(t, span=span, N=nbasis, degree=1),
    name="tent", nbasis=as.integer(nbasis), span=span,
    params=list(N=nbasis, degree=1, span=span)
  )
  class(obj) <- c("Tent_HRF", class(obj))
  obj
}

#' Create Fourier HRF Basis Set
#'
#' Generates an HRF object using Fourier basis functions (sine and cosine pairs)
#' with custom parameters.
#'
#' @param nbasis Number of basis functions (default: 5). Should be even for complete sine-cosine pairs.
#' @param span Temporal window in seconds (default: 24)
#' @return An HRF object of class \code{c("Fourier_HRF", "HRF", "function")}
#' @details 
#' The Fourier basis uses alternating sine and cosine functions with increasing
#' frequencies. This provides a smooth, periodic basis set that can capture
#' oscillatory components in the HRF.
#' @seealso \code{\link{HRF_objects}} for pre-defined HRF objects,
#'   \code{\link{getHRF}} for a unified interface to create HRFs
#' @examples
#' # Create Fourier basis with 8 functions
#' custom_fourier <- hrf_fourier_generator(nbasis = 8)
#' t <- seq(0, 24, by = 0.1)
#' response <- evaluate(custom_fourier, t)
#' matplot(t, response, type = "l", main = "Fourier HRF with 8 basis functions")
#' @export
hrf_fourier_generator <- function(nbasis=5, span=24) {
  obj <- as_hrf(
    f = function(t) hrf_fourier(t, span=span, nbasis=nbasis),
    name="fourier", nbasis=as.integer(nbasis), span=span,
    params=list(nbasis=nbasis, span=span)
  )
  class(obj) <- c("Fourier_HRF", class(obj))
  obj
}

#' Create Daguerre HRF Basis Set
#'
#' Generates an HRF object using Daguerre spherical basis functions with custom parameters.
#' These are orthogonal polynomials that naturally decay to zero.
#'
#' @param nbasis Number of basis functions (default: 3)
#' @param scale Scale parameter for the time axis (default: 4)
#' @return An HRF object of class \code{c("Daguerre_HRF", "HRF", "function")}
#' @details 
#' Daguerre basis functions are orthogonal polynomials on [0,Inf) with respect
#' to the weight function w(x) = x^2 * exp(-x). They are particularly useful
#' for modeling hemodynamic responses as they naturally decay to zero and can
#' capture various response shapes with few parameters.
#' @seealso \code{\link{HRF_objects}} for pre-defined HRF objects,
#'   \code{\link{getHRF}} for a unified interface to create HRFs
#' @examples
#' # Create Daguerre basis with 5 functions
#' custom_dag <- hrf_daguerre_generator(nbasis = 5, scale = 3)
#' t <- seq(0, 24, by = 0.1)
#' response <- evaluate(custom_dag, t)
#' matplot(t, response, type = "l", main = "Daguerre HRF with 5 basis functions")
#' @export
hrf_daguerre_generator <- function(nbasis=3, scale=4) {
  obj <- as_hrf(
    f = function(t) daguerre_basis(t, n_basis=nbasis, scale=scale),
    name="daguerre", nbasis=as.integer(nbasis), span=24,
    params=list(n_basis=nbasis, scale=scale)
  )
  class(obj) <- c("Daguerre_HRF", class(obj))
  obj
}

#' Create FIR HRF Basis Set
#'
#' Generates an HRF object using Finite Impulse Response (FIR) basis functions
#' with custom parameters. Each basis function represents a time bin with a
#' value of 1 in that bin and 0 elsewhere.
#'
#' @param nbasis Number of time bins (default: 12)
#' @param span Temporal window in seconds (default: 24)
#' @return An HRF object of class \code{c("FIR_HRF", "HRF", "function")}
#' @details 
#' The FIR basis divides the time window into \code{nbasis} equal bins.
#' Each basis function is an indicator function for its corresponding bin.
#' This provides maximum flexibility but requires more parameters than
#' smoother basis sets like B-splines.
#' @seealso \code{\link{HRF_objects}} for pre-defined HRF objects,
#'   \code{\link{getHRF}} for a unified interface to create HRFs,
#'   \code{\link{hrf_bspline_generator}} for a smoother alternative
#' @examples
#' # Create FIR basis with 20 bins over 30 seconds
#' custom_fir <- hrf_fir_generator(nbasis = 20, span = 30)
#' t <- seq(0, 30, by = 0.1)
#' response <- evaluate(custom_fir, t)
#' matplot(t, response, type = "l", main = "FIR HRF with 20 time bins")
#' 
#' # Compare to default FIR with 12 bins
#' default_fir <- HRF_FIR
#' response_default <- evaluate(default_fir, t[1:241])  # 24 seconds
#' matplot(t[1:241], response_default, type = "l", 
#'         main = "Default FIR HRF (12 bins over 24s)")
#' @export
hrf_fir_generator <- function(nbasis = 12, span = 24) {
  assertthat::assert_that(
    is.numeric(nbasis) && length(nbasis) == 1 && nbasis >= 1,
    msg = "`nbasis` must be a single positive integer."
  )
  assertthat::assert_that(
    is.numeric(span) && length(span) == 1 && span > 0,
    msg = "`span` must be a single positive number."
  )
  nbasis <- as.integer(nbasis)
  bin_width <- span / nbasis

  f_fir <- function(t) {
    if (!is.numeric(t) || length(t) == 0) {
      return(matrix(0, nrow = 0, ncol = nbasis))
    }
    output_matrix <- matrix(0, nrow = length(t), ncol = nbasis)
    for (i in seq_along(t)) {
      current_t <- t[i]
      if (!is.na(current_t) && current_t >= 0 && current_t < span) {
        bin_index <- if (current_t == 0) 1 else floor(current_t / bin_width) + 1
        bin_index <- min(bin_index, nbasis) 
        output_matrix[i, bin_index] <- 1
      }
    }
    return(output_matrix)
  }

  obj <- as_hrf(
    f = f_fir,
    name = "fir",
    nbasis = nbasis,
    span = span
  )
  # Store params as metadata (the closure already captured these values)
  attr(obj, "params") <- list(nbasis = nbasis, span = span, bin_width = bin_width)
  attr(obj, "param_names") <- c("nbasis", "span", "bin_width")
  class(obj) <- c("FIR_HRF", class(obj))
  obj
}

# Define HRF Registry -----

#' @keywords internal
HRF_REGISTRY <- list(
  spmg1    = HRF_SPMG1,
  spmg2    = HRF_SPMG2,
  spmg3    = HRF_SPMG3,
  gamma    = HRF_GAMMA,
  gaussian = HRF_GAUSSIAN,
  bspline  = hrf_bspline_generator,
  tent     = hrf_tent_generator,
  fourier  = hrf_fourier_generator,
  daguerre = hrf_daguerre_generator,
  fir      = hrf_fir_generator,
  lwu      = hrf_lwu,
  # Aliases
  gam      = HRF_GAMMA,
  bs       = hrf_bspline_generator
)

# getHRF function using the registry (Minimal Version) -----

#' Get HRF by Name
#'
#' Retrieves an HRF by name from the registry and optionally applies decorators.
#' This provides a unified interface for creating both pre-defined HRF objects
#' and custom basis sets with specified parameters.
#'
#' @param name Character string specifying the HRF type. Options include:
#'   \itemize{
#'     \item \code{"spmg1"}, \code{"spmg2"}, \code{"spmg3"} - SPM canonical HRFs
#'     \item \code{"gamma"}, \code{"gaussian"} - Simple parametric HRFs
#'     \item \code{"fir"} - Finite Impulse Response basis
#'     \item \code{"bspline"} or \code{"bs"} - B-spline basis
#'     \item \code{"fourier"} - Fourier basis
#'     \item \code{"daguerre"} - Daguerre spherical basis
#'     \item \code{"tent"} - Tent (linear spline) basis
#'   }
#' @param nbasis Number of basis functions (for basis set types)
#' @param span Temporal window in seconds (default: 24)
#' @param lag Time lag in seconds to apply (default: 0)
#' @param width Block width for block designs (default: 0)
#' @param summate Whether to sum responses in block designs (default: TRUE)
#' @param normalize Whether to normalize the HRF (default: FALSE)
#' @param ... Additional arguments passed to generator functions (e.g., \code{scale} for daguerre)
#' @return An HRF object
#' @details
#' For single HRF types (spmg1, gamma, gaussian), the function returns
#' pre-defined objects. For basis set types (fir, bspline, fourier, daguerre),
#' it calls the appropriate generator function with the specified parameters.
#' @examples
#' # Get pre-defined canonical HRF
#' canonical <- getHRF("spmg1")
#' 
#' # Create custom FIR basis with 20 bins
#' fir20 <- getHRF("fir", nbasis = 20, span = 30)
#' 
#' # Create B-spline basis with lag
#' bs_lag <- getHRF("bspline", nbasis = 8, lag = 2)
#' 
#' # Create blocked Gaussian HRF
#' block_gauss <- getHRF("gaussian", width = 5)
#' @export
getHRF <- function(name = "spmg1", # Default to spmg1
                   nbasis=5, span=24,
                   lag=0, width=0,
                   summate=TRUE, normalize=FALSE, ...) {

  key   <- match.arg(tolower(name), names(HRF_REGISTRY))
  entry <- HRF_REGISTRY[[key]]

  base <- if (inherits(entry, "HRF")) {
            entry # Use pre-defined object
      } else {
            # Call generator, passing nbasis, span, and any relevant ... args
            gen_args <- c(list(nbasis=as.integer(nbasis), span=span), list(...))
            # Only pass args the generator actually accepts
            valid_args <- gen_args[names(gen_args) %in% names(formals(entry))]
            do.call(entry, valid_args)
          }

  # Apply decorators
  if (width != 0) {
      stopifnot(width > 0)
      base <- block_hrf(base, width = width, summate = summate)
  }
  if (lag != 0) {
      base <- lag_hrf(base, lag = lag)
  }
  if (normalize) {
      base <- normalise_hrf(base)
  }

  attr(base, "name") <- key # Set name attribute to the matched registry key
  base
}

#' Evaluate an HRF Object
#'
#' This function evaluates a hemodynamic response function (HRF) object for a given set of time points (grid) and other parameters.
#' It handles both point evaluation (duration=0) and block evaluation (duration > 0).
#'
#' @param x The HRF object (inherits from `HRF` and `function`).
#' @param grid A numeric vector of time points at which to evaluate the HRF.
#' @param amplitude The scaling value for the event (default: 1).
#' @param duration The duration of the event (seconds). If > 0, the HRF is evaluated over this duration (default: 0).
#' @param precision The temporal resolution for evaluating responses when duration > 0 (default: 0.2).
#' @param summate Logical; whether the HRF response should accumulate over the duration (default: TRUE). If FALSE, the maximum response within the duration window is taken (currently only supported for single-basis HRFs).
#' @param normalize Logical; scale output so that the peak absolute value is 1 (default: FALSE). Applied *after* amplitude scaling and duration processing.
#' @param ... Additional arguments (unused).
#' @return A numeric vector or matrix of HRF values at the specified time points.
#' @examples
#' # Evaluate canonical HRF at specific times
#' times <- seq(0, 20, by = 0.5)
#' response <- evaluate(HRF_SPMG1, times)
#' 
#' # Evaluate with amplitude scaling
#' response_scaled <- evaluate(HRF_SPMG1, times, amplitude = 2)
#' 
#' # Evaluate with duration (block design)
#' response_block <- evaluate(HRF_SPMG1, times, duration = 5, summate = TRUE)
#' 
#' # Multi-basis HRF evaluation
#' response_multi <- evaluate(HRF_SPMG3, times)  # Returns 3-column matrix
#' @export
evaluate.HRF <- function(x, grid, amplitude = 1, duration = 0,
                         precision = .2, summate = TRUE, normalize = FALSE, ...) {

  # Validate inputs
  if (!is.numeric(grid) || length(grid) == 0 || anyNA(grid)) {
    stop("`grid` must be a non-empty numeric vector with no NA values.",
         call. = FALSE)
  }
  if (!is.numeric(precision) || length(precision) != 1 || is.na(precision) ||
      precision <= 0) {
    stop("`precision` must be a positive numeric value.", call. = FALSE)
  }

  # Base function incorporating amplitude
  base <- function(g) amplitude * x(g)

  # Evaluate based on duration
  out <- if (duration < precision) {
      # Point evaluation
    base(grid)
  } else {
      # Block evaluation
    offs <- seq(0, duration, by = precision)
      # Evaluate HRF at shifted time points for each offset
      # Use lapply to handle potential matrix output from multi-basis HRFs
      hlist <- lapply(offs, function(o) base(grid - o))
      
      # Check if the result for the first offset is a matrix (multi-basis)
      is_multi_basis <- is.matrix(hlist[[1]])
      
      if (is_multi_basis) {
          # Combine matrices (summation is standard for multi-basis)
          if (summate) {
             Reduce("+", hlist)
    } else {
             # Taking max per-basis-column across offsets is non-standard and complex.
             # Sticking to summation for multi-basis block designs.
             warning("summate=FALSE is not typically used with multi-basis HRFs during block evaluation. Using summation.", call. = FALSE)
             Reduce("+", hlist)
          }
      } else {
          # Single basis HRF: hlist contains vectors, bind them into a matrix
        hmat <- do.call(cbind, hlist)
          if (summate) {
              rowSums(hmat)
          } else {
              # For single basis, take the max across the duration window at each grid point
              apply(hmat, 1, max, na.rm = TRUE) 
              # Alternative: find which offset gives max? apply(hmat, 1, function(vals) vals[which.max(vals)])
          }
      }
  }

  # Apply normalization if requested, handling matrix/vector case
  if (normalize) {
      if (is.matrix(out)) {
          peaks <- apply(out, 2, function(col) max(abs(col), na.rm = TRUE))
          peaks[peaks == 0 | is.na(peaks)] <- 1
          out <- sweep(out, 2, peaks, "/")
      } else {
          peak_val <- max(abs(out), na.rm = TRUE)
          if (!is.na(peak_val) && peak_val != 0) out / peak_val else out
      }
  } else {
    out
  }
}

#' Plot an HRF Object
#'
#' Creates a visualization of an HRF object. For single-basis HRFs, shows the
#' response curve with peak annotation. For multi-basis HRFs (e.g., HRF_SPMG3),
#' shows all basis functions on the same plot.
#'
#' @param x An HRF object
#' @param time Numeric vector of time points. If NULL (default), uses
#'   seq(0, span, by = 0.1) where span is the HRF's span attribute.
#' @param normalize Logical; if TRUE, normalize responses to peak at 1.
#'   Default is FALSE.
#' @param show_peak Logical; if TRUE (default for single-basis HRFs), annotate
#'   the peak time and amplitude on the plot.
#' @param ... Additional arguments passed to underlying plot functions.
#' @return Invisibly returns a data frame with the time and response values
#'   (useful for further customization).
#' @examples
#' # Plot single-basis HRF
#' plot(HRF_SPMG1)
#'
#' # Plot multi-basis HRF
#' plot(HRF_SPMG3)
#'
#' # Plot with normalization
#' plot(HRF_GAMMA, normalize = TRUE)
#'
#' # Custom time range
#' plot(HRF_SPMG1, time = seq(0, 30, by = 0.5))
#' @method plot HRF
#' @export
plot.HRF <- function(x, time = NULL, normalize = FALSE, show_peak = TRUE, ...) {
  span <- attr(x, "span") %||% 24
  if (is.null(time)) {
    time <- seq(0, span, by = 0.1)
  }

  y <- evaluate(x, time, normalize = normalize)
  hrf_name <- attr(x, "name") %||% "HRF"

  if (is.matrix(y)) {
    # Multi-basis HRF
    graphics::matplot(time, y, type = "l", lwd = 1.5, lty = 1,
                      xlab = "Time (s)", ylab = "Response",
                      main = hrf_name, ...)
    nb <- ncol(y)
    graphics::legend("topright", paste("Basis", 1:nb),
                     col = 1:nb, lty = 1, lwd = 1.5, bty = "n")
  } else {
    # Single-basis HRF
    graphics::plot(time, y, type = "l", lwd = 1.5,
                   xlab = "Time (s)", ylab = "Response",
                   main = hrf_name, ...)
    if (show_peak && length(y) > 0) {
      peak_idx <- which.max(y)
      graphics::points(time[peak_idx], y[peak_idx], pch = 19, col = "red")
      graphics::text(time[peak_idx], y[peak_idx],
                     sprintf("Peak: %.1fs", time[peak_idx]),
                     pos = 3, offset = 0.5, col = "red")
    }
  }

  # Return data invisibly for further use
  if (is.matrix(y)) {
    df <- data.frame(time = time, y)
    colnames(df)[-1] <- paste0("basis_", 1:ncol(y))
  } else {
    df <- data.frame(time = time, response = y)
  }
  invisible(df)
}

#' Compare Multiple HRF Functions
#'
#' Creates a comparison plot of multiple HRF objects. This function provides
#' a convenient way to visualize different HRFs on the same plot, with options
#' for normalization and customization. Uses ggplot2 if available for
#' publication-quality figures, otherwise falls back to base R graphics.
#'
#' @param ... HRF objects to compare. Can be passed as individual arguments
#'   or as a named list.
#' @param time Numeric vector of time points. If NULL (default), uses
#'   seq(0, max_span, by = 0.1) where max_span is the maximum span across
#'   all HRFs.
#' @param normalize Logical; if TRUE, normalize all HRFs to peak at 1.
#'   Useful for comparing shapes regardless of amplitude. Default is FALSE.
#' @param labels Character vector of labels for each HRF. If NULL (default),
#'   uses the 'name' attribute of each HRF, or "HRF_1", "HRF_2", etc.
#' @param title Character string for the plot title. If NULL (default),
#'   uses "HRF Comparison".
#' @param subtitle Character string for the plot subtitle. If NULL (default),
#'   no subtitle is shown.
#' @param use_ggplot Logical; if TRUE and ggplot2 is available, use ggplot2
#'   for plotting. If FALSE, use base R graphics. Default is TRUE.
#' @return Invisibly returns a data frame in long format with columns 'time',
#'   'HRF', and 'response'. If use_ggplot is TRUE and ggplot2 is available,
#'   also returns a ggplot object as an attribute 'plot'.
#' @examples
#' # Compare canonical HRFs
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, HRF_GAUSSIAN)
#'
#' # Compare with custom labels
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA,
#'           labels = c("SPM Canonical", "Gamma"))
#'
#' # Normalize for shape comparison
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, HRF_GAUSSIAN,
#'           normalize = TRUE,
#'           title = "HRF Shape Comparison",
#'           subtitle = "All HRFs normalized to peak at 1")
#'
#' # Compare blocked HRFs with different durations
#' hrf_1s <- block_hrf(HRF_SPMG1, width = 1)
#' hrf_3s <- block_hrf(HRF_SPMG1, width = 3)
#' hrf_5s <- block_hrf(HRF_SPMG1, width = 5)
#' plot_hrfs(hrf_1s, hrf_3s, hrf_5s,
#'           labels = c("1s duration", "3s duration", "5s duration"),
#'           title = "Effect of Event Duration on HRF")
#'
#' # Use base R graphics instead of ggplot2
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, use_ggplot = FALSE)
#' @export
plot_hrfs <- function(..., time = NULL, normalize = FALSE, labels = NULL,
                      title = NULL, subtitle = NULL, use_ggplot = TRUE) {
  # Collect HRF objects
  hrfs <- list(...)

  # Handle case where a single list is passed

if (length(hrfs) == 1 && is.list(hrfs[[1]]) && !inherits(hrfs[[1]], "HRF")) {
    hrfs <- hrfs[[1]]
  }

  n_hrfs <- length(hrfs)
  if (n_hrfs == 0) {
    stop("At least one HRF object must be provided")
  }

  # Validate all inputs are HRF objects
  for (i in seq_along(hrfs)) {
    if (!inherits(hrfs[[i]], "HRF")) {
      stop("All arguments must be HRF objects. Argument ", i, " is not an HRF.")
    }
  }

  # Determine time points
  if (is.null(time)) {
    spans <- sapply(hrfs, function(h) attr(h, "span") %||% 24)
    max_span <- max(spans)
    time <- seq(0, max_span, by = 0.1)
  }

  # Generate labels
  if (is.null(labels)) {
    labels <- sapply(hrfs, function(h) attr(h, "name") %||% "HRF")
    # Make unique if needed
    if (any(duplicated(labels))) {
      labels <- paste0(labels, "_", seq_along(labels))
    }
  }
  if (length(labels) != n_hrfs) {
    stop("Length of 'labels' must match number of HRFs")
  }

  # Evaluate all HRFs
  responses <- lapply(hrfs, function(h) {
    y <- evaluate(h, time, normalize = normalize)
    # For multi-basis HRFs, take the first basis (canonical)
    if (is.matrix(y)) {
      y[, 1]
    } else {
      y
    }
  })

  # Build data frame
  df_list <- lapply(seq_along(responses), function(i) {
    data.frame(time = time, HRF = labels[i], response = responses[[i]])
  })
  df <- do.call(rbind, df_list)
  df$HRF <- factor(df$HRF, levels = labels)  # Preserve order

  # Set default title
  if (is.null(title)) {
    title <- "HRF Comparison"
  }

  # Plot
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  if (use_ggplot && has_ggplot) {
    # ggplot2 version - use aes() with unquoted names (standard NSE)
    # Note: using aes_string() is deprecated, but direct aes() with unquoted
    # column names works fine when the columns exist in the data
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = response, color = HRF)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = "Time (seconds)",
        y = if (normalize) "Response (normalized)" else "Response",
        color = "HRF"
      ) +
      ggplot2::theme_minimal()

    print(p)
    attr(df, "plot") <- p
  } else {
    # Base R version
    colors <- grDevices::rainbow(n_hrfs)
    y_range <- range(df$response, na.rm = TRUE)

    graphics::plot(NULL, xlim = range(time), ylim = y_range,
                   xlab = "Time (seconds)",
                   ylab = if (normalize) "Response (normalized)" else "Response",
                   main = title)
    if (!is.null(subtitle)) {
      graphics::mtext(subtitle, side = 3, line = 0.5, cex = 0.8)
    }

    for (i in seq_along(responses)) {
      graphics::lines(time, responses[[i]], col = colors[i], lwd = 1.5)
    }

    graphics::legend("topright", legend = labels, col = colors, lty = 1,
                     lwd = 1.5, bty = "n")
  }

  invisible(df)
}

#' Print an HRF Object
#'
#' Displays a concise summary of an HRF object including its name, number of basis
#' functions, temporal span, and parameters (if any).
#'
#' @param x An HRF object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the HRF object
#' @examples
#' # Print canonical HRF
#' print(HRF_SPMG1)
#'
#' # Print multi-basis HRF
#' print(HRF_SPMG3)
#'
#' # Print Gaussian HRF
#' print(HRF_GAUSSIAN)
#' @method print HRF
#' @export
print.HRF <- function(x, ...) {
  hrf_name <- attr(x, "name") %||% "unnamed"
  nb <- attr(x, "nbasis") %||% 1L
  span <- attr(x, "span") %||% NA
  params <- attr(x, "params")
  param_names <- attr(x, "param_names")

  # Use cat() for reliable output in all contexts (including knitr/vignettes)
  cat("-- HRF:", hrf_name, paste(rep("-", max(1, 50 - nchar(hrf_name))), collapse = ""), "\n")
  cat("   Basis functions:", nb, "\n")
  if (!is.na(span)) {
    cat("   Span:", span, "s\n")
  }

  # Show parameters if they exist
  if (!is.null(params) && length(params) > 0) {
    param_str <- paste(
      names(params),
      sapply(params, function(p) format(p, digits = 4)),
      sep = " = ",
      collapse = ", "
    )
    cat("   Parameters:", param_str, "\n")
  } else if (!is.null(param_names) && length(param_names) > 0) {
    cat("   Parameter names:", paste(param_names, collapse = ", "), "\n")
  }

  invisible(x)
}

# Create pre-defined HRF objects using generators -----

#' @rdname HRF_objects
#' @export
HRF_BSPLINE <- hrf_bspline_generator(nbasis=5, span=24)

#' @rdname HRF_objects  
#' @export
HRF_FIR <- hrf_fir_generator(nbasis=12, span=24)
