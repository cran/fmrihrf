#' @importFrom splines bs
#' @importFrom stats dgamma dnorm quantile
NULL

#' HRF (hemodynamic response function) as a linear function of time
#'
#' The `hrf_time` function computes the value of an HRF, which is a simple linear function of time `t`, when `t` is greater than 0 and less than `maxt`.
#'
#' @param t A numeric value representing time in seconds.
#' @param maxt A numeric value representing the maximum time point in the domain. Default value is 22.
#' @return A numeric value representing the value of the HRF at the given time `t`.
#' @family hrf_functions
#' @export
#' @examples
#' # Compute the HRF value for t = 5 seconds with the default maximum time
#' hrf_val <- hrf_time(5)
#'
#' # Compute the HRF value for t = 5 seconds with a custom maximum time of 30 seconds
#' hrf_val_custom_maxt <- hrf_time(5, maxt = 30)
hrf_time <- function(t, maxt=22) {
  ifelse(t > 0 & t < maxt, t, 0)
}

# hrf_ident
# 
# @param t time in seconds
# @export
hrf_ident <- function(t) {
  ifelse( t == 0, 1, 0)
}

#' B-spline HRF (hemodynamic response function)
#'
#' The `hrf_bspline` function computes the B-spline representation of an HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param span A numeric value representing the temporal window over which the basis set spans. Default value is 20.
#' @param N An integer representing the number of basis functions. Default value is 5.
#' @param degree An integer representing the degree of the spline. Default value is 3.
#' @return A matrix representing the B-spline basis for the HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the B-spline HRF representation for time points from 0 to 20 with 0.5 increments
#' hrfb <- hrf_bspline(seq(0, 20, by = .5), N = 4, degree = 2)
#' @export
#' @importFrom splines bs
#' @param ... Additional arguments passed to `splines::bs`.
hrf_bspline <- function(t, span=24, N=5, degree=3, ...) {
	
	ord <- 1 + degree
	# Check if requested N is sufficient for the degree
	if (N < ord) {
	    warning(paste0("Requested N=", N, " basis functions is less than degree+1=", ord, ". ",
	                   "Using minimum required of ", ord, " basis functions."))
	    # We don't change N here, let splines::bs handle the df inconsistency if needed,
	    # but the warning informs the user.
	}
	
	nIknots <- N - ord + 1
	if (nIknots < 0) {
		nIknots <- 0
		#warning("'df' was too small; have used  ", ord - (1 - intercept))
	}
	
	knots <- if (nIknots > 0) {
				knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
				stats::quantile(seq(0,span), knots)
			} else {
				0
			}
	
	t <- as.numeric(t)
	in_support <- !is.na(t) & t >= 0 & t <= span
	t_eval <- t
	t_eval[!in_support] <- 0

	basis <- splines::bs(t_eval, df=N, knots=knots, degree=degree, Boundary.knots=c(0,span),...)
	if (any(!in_support)) {
		basis[!in_support, ] <- 0
	}
	basis
}


#' Gamma HRF (hemodynamic response function)
#'
#' The `hrf_gamma` function computes the gamma density-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param shape A numeric value representing the shape parameter for the gamma probability density function. Default value is 6.
#' @param rate A numeric value representing the rate parameter for the gamma probability density function. Default value is 1.
#' @return A numeric vector representing the gamma HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the gamma HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_gamma_vals <- hrf_gamma(seq(0, 20, by = .5), shape = 6, rate = 1)
#' @export
hrf_gamma <- function(t, shape=6, rate=1) {
  stats::dgamma(t, shape=shape, rate=rate)
}

#' Gaussian HRF (hemodynamic response function)
#'
#' The `hrf_gaussian` function computes the Gaussian density-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param mean A numeric value representing the mean of the Gaussian probability density function. Default value is 6.
#' @param sd A numeric value representing the standard deviation of the Gaussian probability density function. Default value is 2.
#' @return A numeric vector representing the Gaussian HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the Gaussian HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_gaussian_vals <- hrf_gaussian(seq(0, 20, by = .5), mean = 6, sd = 2)
#' @export
hrf_gaussian <- function(t, mean=6, sd=2) {
	stats::dnorm(t, mean=mean, sd=sd)
}



#' Mexican Hat HRF (hemodynamic response function)
#'
#' The `hrf_mexhat` function computes the Mexican hat wavelet-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param mean A numeric value representing the mean of the Mexican hat wavelet. Default value is 6.
#' @param sd A numeric value representing the standard deviation of the Mexican hat wavelet. Default value is 2.
#' @return A numeric vector representing the Mexican hat wavelet-based HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the Mexican hat HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_mexhat_vals <- hrf_mexhat(seq(0, 20, by = .5), mean = 6, sd = 2)
#' @export
hrf_mexhat <- function(t, mean = 6, sd = 2) {
  t0 <- t - mean
  a <- (1 - (t0 / sd)^2) * exp(-t0^2 / (2 * sd^2))
  scale <- sqrt(2 / (3 * sd * pi^(1/4)))
  return(scale * a)
}

#' hrf_spmg1
#'
#' A hemodynamic response function based on the SPM canonical double gamma parameterization.
#'
#' This function models the hemodynamic response using the canonical double gamma parameterization
#' in the SPM software. The HRF is defined by a linear combination of two gamma functions with
#' different exponents (P1 and P2) and amplitudes (A1 and A2). It is commonly used in fMRI data
#' analysis to estimate the BOLD (blood-oxygen-level-dependent) signal changes associated with
#' neural activity.
#'
#' @param t A vector of time points.
#' @param P1 The first exponent parameter (default: 5).
#' @param P2 The second exponent parameter (default: 15).
#' @param A1 Amplitude scaling factor for the positive gamma function component; normally fixed at .0833
#' @return A vector of HRF values at the given time points.
#' @family hrf_functions
#' @export
#' @examples
#' # Generate a time vector
#' time_points <- seq(0, 30, by=0.1)
#' # Compute the HRF values using the SPM canonical double gamma parameterization
#' hrf_values <- hrf_spmg1(time_points)
#' # Plot the HRF values
#' plot(time_points, hrf_values, type='l', main='SPM Canonical Double Gamma HRF')
hrf_spmg1 <- function(t, P1=5, P2=15,A1=.0833) {
 	ifelse(t < 0, 0, exp(-t)*(A1*t^P1 - 1.274527e-13*t^P2))
	
}


# Fast analytic first derivative for hrf_spmg1
#' @keywords internal
#' @noRd
hrf_spmg1_deriv <- function(t, P1 = 5, P2 = 15, A1 = .0833) {
  C <- 1.274527e-13
  ret <- numeric(length(t))
  pos <- t >= 0
  if (any(pos)) {
    t_pos <- t[pos]
    ret[pos] <- exp(-t_pos) * (A1 * t_pos^(P1 - 1) * (P1 - t_pos) -
                                 C   * t_pos^(P2 - 1) * (P2 - t_pos))
  }
  ret
}

# Fast analytic second derivative for hrf_spmg1
#' @keywords internal
#' @noRd
hrf_spmg1_second_deriv <- function(t, P1 = 5, P2 = 15, A1 = .0833) {
  C <- 1.274527e-13
  ret <- numeric(length(t))
  pos <- t >= 0
  if (any(pos)) {
    t_pos <- t[pos]
    # Let D1 = A1 * t^(P1-1) * (P1 - t) and D2 = C * t^(P2-1) * (P2 - t)
    D1 <- A1 * t_pos^(P1 - 1) * (P1 - t_pos)
    D2 <- C   * t_pos^(P2 - 1) * (P2 - t_pos)
    # Their derivatives:
    D1_prime <- A1 * ((P1 - 1) * t_pos^(P1 - 2) * (P1 - t_pos) - t_pos^(P1 - 1))
    D2_prime <- C   * ((P2 - 1) * t_pos^(P2 - 2) * (P2 - t_pos) - t_pos^(P2 - 1))
    ret[pos] <- exp(-t_pos) * (D1_prime - D2_prime - (D1 - D2))
  }
  ret
}

#' hrf_sine
#'
#' A hemodynamic response function using the Sine Basis Set.
#'
#' @param t A vector of times.
#' @param span The temporal window over which the basis sets span (default: 24).
#' @param N The number of basis functions (default: 5).
#' @return A matrix of sine basis functions.
#' @family hrf_functions
#' @export
#' @examples
#' hrf_sine_basis <- hrf_sine(seq(0, 20, by = 0.5), N = 4)
hrf_sine <- function(t, span = 24, N = 5) {
  t <- as.numeric(t)
  in_support <- !is.na(t) & t >= 0 & t <= span
  sine_basis <- vapply(1:N, function(n) {
    sin(2 * pi * n * t / span)
  }, numeric(length(t)))
  if (any(!in_support)) {
    sine_basis[!in_support, ] <- 0
  }
  return(sine_basis)
}

#' hrf_inv_logit
#'
#' A hemodynamic response function using the difference of two Inverse Logit functions.
#'
#' @param t A vector of times.
#' @param mu1 The time-to-peak for the rising phase (mean of the first logistic function).
#' @param s1 The width (slope) of the first logistic function.
#' @param mu2 The time-to-peak for the falling phase (mean of the second logistic function).
#' @param s2 The width (slope) of the second logistic function.
#' @param lag The time delay (default: 0).
#' @return A vector of the difference of two Inverse Logit HRF values.
#' @family hrf_functions
#' @export
#' @examples
#' hrf_inv_logit_basis <- hrf_inv_logit(seq(0, 20, by = 0.5), mu1 = 6, s1 = 1, mu2 = 16, s2 = 1)
hrf_inv_logit <- function(t, mu1 = 6, s1 = 1, mu2 = 16, s2 = 1, lag = 0) {
  inv_logit1 <- 1 / (1 + exp(-(t - lag - mu1) / s1))
  inv_logit2 <- 1 / (1 + exp(-(t - lag - mu2) / s2))
  return(inv_logit1 - inv_logit2)
}


#' Hemodynamic Response Function with Half-Cosine Basis
#'
#' This function models a hemodynamic response function (HRF) using four half-period cosine basis functions.
#' The HRF consists of an initial dip, a rise to peak, a fall and undershoot, and a recovery to the baseline.
#'
#' @param t A vector of time values.
#' @param h1 Duration of the initial dip in seconds.
#' @param h2 Duration of the rise to peak in seconds.
#' @param h3 Duration of the fall and undershoot in seconds.
#' @param h4 Duration of the recovery to baseline in seconds.
#' @param f1 Height of the starting point.
#' @param f2 Height of the end point.
#' @return A vector of HRF values corresponding to the input time values.
#' @references Woolrich, M. W., Behrens, T. E., & Smith, S. M. (2004). Constrained linear basis sets for HRF modelling using Variational Bayes. NeuroImage, 21(4), 1748-1761.
#'
#' Half-cosine HRF
#'
#' Creates a hemodynamic response function using half-cosine segments.
#' The function consists of four phases controlled by h1-h4 parameters,
#' with transitions between baseline (f1) and peak (1) and final (f2) levels.
#'
#' @param t Time points at which to evaluate the HRF
#' @param h1 Duration of initial fall from f1 to 0 (default: 1)
#' @param h2 Duration of rise from 0 to 1 (default: 5)
#' @param h3 Duration of fall from 1 to 0 (default: 7)
#' @param h4 Duration of final rise from 0 to f2 (default: 7)
#' @param f1 Initial baseline level (default: 0)
#' @param f2 Final baseline level (default: 0)
#' @return Numeric vector of HRF values at time points t
#' @examples
#' # Standard half-cosine HRF
#' t <- seq(0, 30, by = 0.1)
#' hrf <- hrf_half_cosine(t)
#' plot(t, hrf, type = "l", main = "Half-cosine HRF")
#' 
#' # Modified shape with undershoot
#' hrf_under <- hrf_half_cosine(t, h1 = 1, h2 = 4, h3 = 6, h4 = 8, f2 = -0.2)
#' lines(t, hrf_under, col = "red")
#' @keywords internal
#' @noRd
.hrf_half_cosine <- function(t, h1=1, h2=5, h3=7,h4=7, f1=0, f2=0) {
  rising_half_cosine <- function(t, f1, t0, w) {
    return(f1/2 * (1 - cos(pi * (t - t0) / w)))
  }
  
  falling_half_cosine <- function(t, f1, t0, w) {
    return(f1/2 * (1 + cos(pi * (t - t0) / w)))
  }
  
  ret = dplyr::case_when(
    t < 0 ~ 0,
    t <= h1 ~ falling_half_cosine(t, f1, 0, h1),
    (t > h1) & t <= (h1+h2) ~ rising_half_cosine(t, 1, h1, h2),
    (t > (h1+h2)) & t <= (h1+h2+h3) ~ falling_half_cosine(t,1,(h1+h2), h3),
    (t > (h1+h2+h3)) & t <= (h1+h2+h3+h4) ~ rising_half_cosine(t,f2,(h1+h2+h3), h4),
    (t > h1+h2+h3+h4) ~ f2,
  )
  return(ret)
}

#' Half-cosine HRF 
#'
#' Segments: 0->f1 (h1), f1->1 (h2), 1->f2 (h3), f2->0 (h4).
#' Negative f1 gives an initial dip; negative f2 gives an undershoot.
#' Peak is at t = h1 + h2 (amplitude 1 by construction).
#'
#' @param t Numeric vector of times (s)
#' @param h1,h2,h3,h4 Segment durations (s). Must be > 0.
#' @param f1 Initial dip level (default 0), typically in [-0.2, 0]
#' @param f2 Undershoot level (default 0), typically in [-0.3, 0]
#' @return Numeric vector same length as t
#' @examples
#' t <- seq(0, 30, by = 0.1)
#' y <- hrf_half_cosine(t)
#' @export
hrf_half_cosine <- function(t, h1=1, h2=5, h3=7, h4=7, f1=0, f2=0) {
  stopifnot(h1 > 0, h2 > 0, h3 > 0, h4 > 0)
  trans <- function(tt, a, b, t0, w) a + 0.5*(b - a) * (1 - cos(pi * (tt - t0)/w))
  t1 <- h1; t2 <- h1 + h2; t3 <- h1 + h2 + h3; t4 <- t3 + h4

  out <- numeric(length(t))
  # segment 1: 0 -> f1
  idx <- t >= 0 & t <= t1
  out[idx] <- trans(t[idx], 0,  f1, 0,   h1)
  # segment 2: f1 -> 1
  idx <- t >  t1 & t <= t2
  out[idx] <- trans(t[idx], f1, 1,  t1,  h2)
  # segment 3: 1 -> f2 (undershoot)
  idx <- t >  t2 & t <= t3
  out[idx] <- trans(t[idx], 1,  f2, t2,  h3)
  # segment 4: f2 -> 0 (recovery)
  idx <- t >  t3 & t <= t4
  out[idx] <- trans(t[idx], f2, 0,  t3,  h4)
  # before/after window: 0
  out[t < 0 | t > t4] <- 0
  out
}

#' Fourier basis for HRF modeling
#'
#' Generates a set of Fourier basis functions (sine and cosine pairs) over a given span.
#'
#' @param t A vector of time points.
#' @param span The temporal window over which the basis functions span (default: 24).
#' @param nbasis The number of basis functions (default: 5). Should be even for full sine-cosine pairs.
#' @return A matrix of Fourier basis functions with nbasis columns.
#' @examples
#' # Create Fourier basis with 5 functions
#' t <- seq(0, 24, by = 0.5)
#' basis <- hrf_fourier(t, span = 24, nbasis = 5)
#' matplot(t, basis, type = "l", main = "Fourier Basis Functions")
#' @export
hrf_fourier <- function(t, span = 24, nbasis = 5) {
  t <- as.numeric(t)
  in_support <- !is.na(t) & t >= 0 & t <= span
  freqs <- ceiling(seq_len(nbasis) / 2)
  basis <- vapply(seq_len(nbasis), function(k) {
    n <- freqs[k]
    if (k %% 2 == 1) {
      sin(2 * pi * n * t / span)
    } else {
      cos(2 * pi * n * t / span)
    }
  }, numeric(length(t)))
  if (any(!in_support)) {
    basis[!in_support, ] <- 0
  }
  return(basis)
}



#' HRF Toeplitz Matrix
#' 
#' @description
#' Create a Toeplitz matrix for hemodynamic response function (HRF) convolution.
#' 
#' @param hrf The hemodynamic response function.
#' @param time A numeric vector representing the time points.
#' @param len The length of the output Toeplitz matrix.
#' @param sparse Logical, if TRUE, the output Toeplitz matrix is returned as a sparse matrix (default: FALSE).
#' 
#' @return A Toeplitz matrix for HRF convolution.
#' @examples
#' # Create HRF and time points
#' hrf_fun <- function(t) hrf_spmg1(t)
#' times <- seq(0, 30, by = 1)
#' 
#' # Create Toeplitz matrix
#' H <- hrf_toeplitz(hrf_fun, times, len = 50)
#' 
#' # Create sparse version
#' H_sparse <- hrf_toeplitz(hrf_fun, times, len = 50, sparse = TRUE)
#' @export
hrf_toeplitz <- function(hrf, time, len, sparse=FALSE) {
  hreg <- hrf(time)
  padding <- len - length(hreg)
  H <- pracma::Toeplitz(c(hreg, rep(0, padding)), c(hreg[1], rep(0, len-1)))
  H <- Matrix::Matrix(H, sparse=sparse)
  H
}


#' Generate Daguerre spherical basis functions
#' 
#' @description
#' Creates a set of Daguerre spherical basis functions. These are orthogonal 
#' polynomials on [0,Inf) with respect to the weight function w(x) = x^2 * exp(-x).
#' They are particularly useful for modeling hemodynamic responses as they naturally
#' decay to zero and can capture various response shapes.
#'
#' @param t Time points at which to evaluate the basis functions
#' @param n_basis Number of basis functions to generate (default: 3)
#' @param scale Scale parameter for the time axis (default: 1)
#' @return A matrix with columns containing the basis functions
#' @keywords internal
#' @noRd
daguerre_basis <- function(t, n_basis = 3, scale = 1) {
  # Scale time
  x <- t/scale
  
  # Initialize basis matrix
  basis <- matrix(0, length(x), n_basis)
  
  # First basis function (n=0)
  basis[,1] <- exp(-x/2)
  
  if(n_basis > 1) {
    # Second basis function (n=1)
    basis[,2] <- (1 - x) * exp(-x/2)
  }
  
  if(n_basis > 2) {
    # Higher order basis functions using recurrence relation
    for(n in 3:n_basis) {
      k <- n - 1
      basis[,n] <- ((2*k - 1 - x) * basis[,n-1] - (k - 1) * basis[,n-2]) / k
    }
  }
  
  # Normalize basis functions
  for(i in 1:n_basis) {
    # Avoid division by zero if a basis function is all zero
    max_abs_val <- max(abs(basis[,i]))
    if (max_abs_val > 1e-10) {
      basis[,i] <- basis[,i] / max_abs_val
    }
  }
  
  basis
}

#' Lag-Width-Undershoot (LWU) HRF
#'
#' Computes the Lag-Width-Undershoot (LWU) hemodynamic response function.
#' This model uses two Gaussian components to model the main response and an optional undershoot.
#' 
#' The LWU model formula combines a positive Gaussian peak and a negative undershoot:
#' h(t; tau, sigma, rho) = exp(-(t-tau)^2/(2*sigma^2)) - rho * exp(-(t-tau-2*sigma)^2/(2*(1.6*sigma)^2))
#'
#' @param t A numeric vector of time points (in seconds).
#' @param tau Lag of the main Gaussian component (time-to-peak of the positive lobe, in seconds). Default: 6.
#' @param sigma Width (standard deviation) of the main Gaussian component (in seconds). Must be > 0.05. Default: 2.5.
#' @param rho Amplitude of the undershoot Gaussian component, relative to the main component. Must be between 0 and 1.5. Default: 0.35.
#' @param normalize Character string specifying normalization type. Either "none" for no normalization (default) or "height" to scale the HRF so its maximum absolute value is 1.
#' @return A numeric vector representing the LWU HRF values at the given time points `t`.
#' @family hrf_functions
#' @export
#' @examples
#' t_points <- seq(0, 30, by = 0.1)
#'
#' # Default LWU HRF
#' lwu_default <- hrf_lwu(t_points)
#' plot(t_points, lwu_default, type = "l", main = "LWU HRF (Default Params)", ylab = "Amplitude")
#'
#' # LWU HRF with no undershoot
#' lwu_no_undershoot <- hrf_lwu(t_points, rho = 0)
#' lines(t_points, lwu_no_undershoot, col = "blue")
#'
#' # LWU HRF with a wider main peak and larger undershoot
#' lwu_custom <- hrf_lwu(t_points, tau = 7, sigma = 1.5, rho = 0.5)
#' lines(t_points, lwu_custom, col = "red")
#' legend("topright", c("Default", "No Undershoot (rho=0)", "Custom (tau=7, sigma=1.5, rho=0.5)"),
#'        col = c("black", "blue", "red"), lty = 1, cex = 0.8)
#'
#' # Height-normalized HRF
#' lwu_normalized <- hrf_lwu(t_points, tau = 6, sigma = 1, rho = 0.35, normalize = "height")
#' plot(t_points, lwu_normalized, type = "l", main = "Height-Normalized LWU HRF", ylab = "Amplitude")
#' abline(h = c(-1, 1), lty = 2, col = "grey") # Max absolute value should be 1
hrf_lwu <- function(t, tau = 6, sigma = 2.5, rho = 0.35, normalize = "none") {
  assertthat::assert_that(is.numeric(t), msg = "`t` must be numeric.")
  assertthat::assert_that(is.numeric(tau) && length(tau) == 1, msg = "`tau` must be a single numeric value.")
  assertthat::assert_that(is.numeric(sigma) && length(sigma) == 1, msg = "`sigma` must be a single numeric value.")
  assertthat::assert_that(sigma > 0.05, msg = "`sigma` must be > 0.05.")
  assertthat::assert_that(is.numeric(rho) && length(rho) == 1, msg = "`rho` must be a single numeric value.")
  assertthat::assert_that(rho >= 0 && rho <= 1.5, msg = "`rho` must be between 0 and 1.5.")
  assertthat::assert_that(normalize %in% c("none", "height", "area"),
                        msg = "`normalize` must be one of 'none', 'height', or 'area'.")

  if (normalize == "area") {
    warning("`normalize = \"area\"` is not yet fully implemented for hrf_lwu and will behave like `normalize = \"none\"`. Area normalization typically requires numerical integration and careful definition of the integration window for HRFs.", call. = FALSE)
    normalize <- "none"
  }

  # Main positive Gaussian component
  term1_exponent <- -((t - tau)^2) / (2 * sigma^2)
  term1 <- exp(term1_exponent)

  # Undershoot Gaussian component
  # tau_u = tau + 2*sigma (peak of undershoot relative to stimulus onset)
  # sigma_u = 1.6*sigma (width of undershoot)
  term2_exponent <- -((t - (tau + 2 * sigma))^2) / (2 * (1.6 * sigma)^2)
  term2 <- rho * exp(term2_exponent)

  response <- term1 - term2

  if (normalize == "height") {
    max_abs_val <- max(abs(response), na.rm = TRUE)
    if (max_abs_val > 1e-10) { # Avoid division by zero or tiny numbers
      response <- response / max_abs_val
    }
  }

  return(response)
}


#' Boxcar HRF (No Hemodynamic Delay)
#'
#' Creates a simple boxcar (step function) HRF that is constant within a time
#' window starting at t=0 and zero outside. Unlike traditional HRFs, this has
#' no hemodynamic delay - it represents an instantaneous response.
#'
#' When used in a GLM, the estimated coefficient represents a (weighted) average
#' of the data within the specified time window. If \code{normalize = TRUE}, the
#' coefficient directly estimates the mean signal in that window.
#'
#' For delayed windows (not starting at t=0), use \code{\link{lag_hrf}} to shift
#' the boxcar in time.
#'
#' @section Note on durations:
#' The \code{width} is fixed when the HRF is created. The \code{duration}
#' parameter in \code{\link{regressor}()} does \strong{not} modify the boxcar
#' width---it controls how long the neural input is sustained (which then gets
#' convolved with this HRF). For trial-varying boxcar widths, use a list of HRFs:
#' \preformatted{
#' widths <- c(4, 6, 8)
#' hrfs <- lapply(widths, function(w) hrf_boxcar(width = w, normalize = TRUE))
#' reg <- regressor(onsets = c(0, 20, 40), hrf = hrfs)
#' }
#'
#' @param width Duration of the boxcar window in seconds.
#' @param amplitude Height of the boxcar (default: 1).
#' @param normalize Logical; if \code{TRUE}, the boxcar is scaled so that its
#'   integral equals 1 (i.e., amplitude = 1/width). This makes the regression
#'   coefficient interpretable as the mean signal in the window.
#'   Default is \code{FALSE}.
#' @return An HRF object that can be used with \code{regressor()} and other
#'   fmrihrf functions.
#' @family hrf_functions
#' @seealso \code{\link{hrf_weighted}} for weighted/shaped boxcars,
#'   \code{\link{lag_hrf}} to shift the window in time
#' @export
#' @examples
#' # Simple boxcar of 5 seconds width
#' hrf1 <- hrf_boxcar(width = 5)
#' t <- seq(-1, 10, by = 0.1)
#' plot(t, evaluate(hrf1, t), type = "s", main = "Simple Boxcar HRF")
#'
#' # Normalized boxcar - coefficient will estimate mean signal in window
#' hrf2 <- hrf_boxcar(width = 5, normalize = TRUE)
#' # integral is now 1, so beta estimates mean(Y[0:5])
#'
#' # Use in a regressor with trial-varying widths
#' hrf_short <- hrf_boxcar(width = 4, normalize = TRUE)
#' hrf_long <- hrf_boxcar(width = 8, normalize = TRUE)
#' reg <- regressor(onsets = c(0, 20), hrf = list(hrf_short, hrf_long))
#'
#' # For delayed windows, use lag_hrf decorator
#' hrf_delayed <- lag_hrf(hrf_boxcar(width = 5), lag = 10)  # Window from 10-15s
hrf_boxcar <- function(width, amplitude = 1, normalize = FALSE) {
  assertthat::assert_that(
    is.numeric(width) && length(width) == 1 && width > 0,
    msg = "`width` must be a single positive numeric value."
  )
  assertthat::assert_that(
    is.numeric(amplitude) && length(amplitude) == 1,
    msg = "`amplitude` must be a single numeric value."
  )

  if (normalize) {
    amplitude <- 1 / width
  }

  f <- function(t) {
    ifelse(t >= 0 & t < width, amplitude, 0)
  }

  as_hrf(f,
         name = sprintf("boxcar[%.2g]", width),
         nbasis = 1L,
         span = width,
         params = list(width = width, amplitude = amplitude, normalize = normalize))
}


#' Weighted HRF (No Hemodynamic Delay)
#'
#' Creates a flexible weighted HRF starting at t=0 with user-specified weights.
#' Unlike traditional HRFs, this has no built-in hemodynamic delay - it directly
#' maps weights to time points, allowing for arbitrary temporal response shapes.
#'
#' This is useful for extracting weighted averages of data at specific time points.
#' When \code{normalize = TRUE} and the HRF is used in a GLM, the estimated
#' coefficient represents a weighted mean of the data at the specified times.
#'
#' There are two ways to specify the temporal structure:
#' \enumerate{
#'   \item \code{width + weights}: Weights are evenly spaced from 0 to \code{width}
#'   \item \code{times + weights}: Explicit time points for each weight (relative to t=0)
#' }
#'
#' For delayed windows (not starting at t=0), use \code{\link{lag_hrf}} to shift
#' the weighted HRF in time.
#'
#' @section Note on durations:
#' The temporal structure (\code{width} or \code{times}) is fixed when the HRF
#' is created. The \code{duration} parameter in \code{\link{regressor}()} does
#' \strong{not} modify the weighted HRF's structure---it controls how long the
#' neural input is sustained (which then gets convolved with this HRF). For
#' trial-varying weighted HRFs, use a list of HRFs:
#' \preformatted{
#' hrf_early <- hrf_weighted(width = 6, weights = c(1, 1, 0, 0), normalize = TRUE)
#' hrf_late <- hrf_weighted(width = 6, weights = c(0, 0, 1, 1), normalize = TRUE)
#' reg <- regressor(onsets = c(0, 20), hrf = list(hrf_early, hrf_late))
#' }
#'
#' @param weights Numeric vector of weights. Required.
#' @param width Total duration of the window in seconds. If provided without
#'   \code{times}, weights are evenly spaced from 0 to \code{width}.
#' @param times Numeric vector of time points (in seconds, relative to t=0) where
#'   weights are specified. Must be strictly increasing and start at 0 for
#'   consistency with other HRFs. If provided, \code{width} is ignored.
#' @param method Interpolation method between time points:
#'   \describe{
#'     \item{"constant"}{Step function - weight is constant until the next time
#'       point (default). Good for discrete time bins.
#'     }
#'     \item{"linear"}{Linear interpolation between points. Good for smooth
#'       weight transitions.
#'     }
#'   }
#' @param normalize Logical; if \code{TRUE}, weights are scaled so they sum to 1
#'   (for \code{method = "constant"}) or integrate to 1 (for \code{method = "linear"}).
#'   This makes the regression coefficient interpretable as a weighted mean.
#'   Default is \code{FALSE}.
#' @return An HRF object that can be used with \code{regressor()} and other
#'   fmrihrf functions.
#' @family hrf_functions
#' @seealso \code{\link{hrf_boxcar}} for simple uniform boxcars,
#'   \code{\link{lag_hrf}} to shift the window in time,
#'   \code{\link{empirical_hrf}} for HRFs from measured data
#' @export
#' @examples
#' # Simple: 6s window with 4 evenly-spaced weights (at 0, 2, 4, 6s)
#' hrf1 <- hrf_weighted(width = 6, weights = c(0.2, 0.5, 0.8, 0.3))
#' t <- seq(-1, 10, by = 0.1)
#' plot(t, evaluate(hrf1, t), type = "s", main = "Weighted HRF (width + weights)")
#'
#' # Explicit times for precise control
#' hrf2 <- hrf_weighted(
#'   times = c(0, 1, 3, 5, 6),
#'   weights = c(0.1, 0.5, 0.8, 0.5, 0.1),
#'   method = "linear"
#' )
#' plot(t, evaluate(hrf2, t), type = "l", main = "Smooth Weighted HRF")
#'
#' # Normalized weights - coefficient estimates weighted mean of signal
#' hrf3 <- hrf_weighted(
#'   width = 8,
#'   weights = c(1, 2, 2, 1),
#'   normalize = TRUE
#' )
#'
#' # Trial-varying weighted HRFs
#' hrf_early <- hrf_weighted(width = 6, weights = c(1, 1, 0, 0), normalize = TRUE)
#' hrf_late <- hrf_weighted(width = 6, weights = c(0, 0, 1, 1), normalize = TRUE)
#' reg <- regressor(onsets = c(0, 20), hrf = list(hrf_early, hrf_late))
#'
#' # For delayed windows, use lag_hrf
#' hrf_delayed <- lag_hrf(hrf_weighted(width = 5, weights = c(1, 2, 1)), lag = 10)
hrf_weighted <- function(weights, width = NULL, times = NULL,
                         method = c("constant", "linear"), normalize = FALSE) {
  method <- match.arg(method)

  assertthat::assert_that(
    is.numeric(weights) && length(weights) >= 2,
    msg = "`weights` must be a numeric vector with at least 2 elements."
  )

  # Determine times: explicit times take precedence, otherwise generate from width
 if (!is.null(times)) {
    assertthat::assert_that(
      is.numeric(times) && length(times) >= 2,
      msg = "`times` must be a numeric vector with at least 2 elements."
    )
    assertthat::assert_that(
      length(times) == length(weights),
      msg = "`times` and `weights` must have the same length."
    )
    assertthat::assert_that(
      all(diff(times) > 0),
      msg = "`times` must be strictly increasing."
    )
    assertthat::assert_that(
      times[1] >= 0,
      msg = "`times` must start at 0 or later (HRFs are relative to event onset)."
    )
  } else if (!is.null(width)) {
    assertthat::assert_that(
      is.numeric(width) && length(width) == 1 && width > 0,
      msg = "`width` must be a single positive numeric value."
    )
    # Generate evenly spaced times from 0 to width
    n_weights <- length(weights)
    times <- seq(0, width, length.out = n_weights)
  } else {
    stop("Either `width` or `times` must be provided.")
  }

  # Normalize weights if requested
  if (normalize) {
    if (method == "constant") {
      # For step function, normalize so weights sum to 1
      # (each weight applies to its interval)
      weight_sum <- sum(weights[-length(weights)])  # last weight has no interval
      if (abs(weight_sum) > 1e-10) {
        # Scale all weights proportionally
        weights <- weights / weight_sum
      }
    } else {
      # For linear interpolation, normalize so integral = 1
      # Use trapezoidal rule to estimate integral
      intervals <- diff(times)
      avg_weights <- (weights[-length(weights)] + weights[-1]) / 2
      integral <- sum(intervals * avg_weights)
      if (abs(integral) > 1e-10) {
        weights <- weights / integral
      }
    }
  }

  # Create the interpolation function
  if (method == "linear") {
    f <- stats::approxfun(times, weights, yleft = 0, yright = 0, rule = 1)
  } else {
    # Piecewise constant (step function)
    # For step function, each weight applies from times[i] to times[i+1]
    f <- stats::approxfun(times, weights, yleft = 0, yright = 0,
                          method = "constant", rule = 1)
  }

  # Name reflects how it was specified
  hrf_name <- if (!is.null(width) && is.null(times)) {
    sprintf("weighted[w=%.2g, %d wts]", max(times), length(weights))
  } else {
    sprintf("weighted[%d pts, %s]", length(times), method)
  }

  as_hrf(f,
         name = hrf_name,
         nbasis = 1L,
         span = max(times),
         params = list(times = times, weights = weights, width = width,
                       method = method, normalize = normalize))
}


#' LWU HRF Basis for Taylor Expansion
#'
#' Constructs the basis set for the Lag-Width-Undershoot (LWU) HRF model,
#' intended for Taylor expansion-based fitting. The basis consists of the
#' LWU HRF evaluated at a given expansion point \code{theta0}, and its
#' partial derivatives with respect to its parameters (tau, sigma, rho).
#'
#' @param theta0 A numeric vector of length 3 specifying the expansion point
#'   \code{c(tau0, sigma0, rho0)} for the LWU parameters.
#' @param t A numeric vector of time points (in seconds) at which to evaluate the basis.
#' @param normalize_primary Character string, one of \code{"none"} or \code{"height"}.
#'   If \code{"height"}, the primary HRF column (\code{h0(t)}) is normalized to have a
#'   peak absolute value of 1. For Taylor expansion fitting as described in Fit_LRU.md,
#'   this should typically be \code{"none"} as the scaling is absorbed by the beta coefficient.
#'   Default is \code{"none"}.
#' @return A numeric matrix of dimension \code{length(t) x 4}. The columns represent:
#'   \itemize{
#'     \item \code{h0}: LWU HRF evaluated at theta0
#'     \item \code{d_tau}: Partial derivative with respect to tau at theta0
#'     \item \code{d_sigma}: Partial derivative with respect to sigma at theta0
#'     \item \code{d_rho}: Partial derivative with respect to rho at theta0
#'   }
#' @family hrf_functions
#' @seealso \code{\link{hrf_lwu}}, \code{\link[numDeriv]{grad}}
#' @export
#' @importFrom numDeriv grad
#' @examples
#' t_points <- seq(0, 30, by = 0.5)
#' theta0_default <- c(tau = 6, sigma = 1, rho = 0.35)
#'
#' # Generate the basis set
#' lwu_basis <- hrf_basis_lwu(theta0_default, t_points)
#' dim(lwu_basis) # Should be length(t_points) x 4
#' head(lwu_basis)
#'
#' # Plot the basis functions
#' matplot(t_points, lwu_basis, type = "l", lty = 1,
#'         main = "LWU HRF Basis Functions", ylab = "Value", xlab = "Time (s)")
#' legend("topright", colnames(lwu_basis), col = 1:4, lty = 1, cex = 0.8)
#'
#' # Example with primary HRF normalization (not typical for Taylor fitting step)
#' lwu_basis_norm_h0 <- hrf_basis_lwu(theta0_default, t_points, normalize_primary = "height")
#' plot(t_points, lwu_basis_norm_h0[,1], type="l", main="Normalized h0 in Basis")
#' max(abs(lwu_basis_norm_h0[,1])) # Should be 1
hrf_basis_lwu <- function(theta0, t, normalize_primary = "none") {
  assertthat::assert_that(is.numeric(theta0) && length(theta0) == 3,
                        msg = "`theta0` must be a numeric vector of length 3: c(tau, sigma, rho).")
  names(theta0) <- c("tau", "sigma", "rho") # Ensure names for numDeriv::grad
  assertthat::assert_that(is.numeric(t), msg = "`t` must be numeric.")
  assertthat::assert_that(normalize_primary %in% c("none", "height"),
                        msg = "`normalize_primary` must be one of 'none' or 'height'.")

  # Safety checks for sigma0 and rho0 from theta0, consistent with hrf_lwu
  assertthat::assert_that(theta0["sigma"] > 0.05, msg = "sigma in `theta0` must be > 0.05.")
  assertthat::assert_that(theta0["rho"] >= 0 && theta0["rho"] <= 1.5,
                        msg = "rho in `theta0` must be between 0 and 1.5.")

  # Function to pass to numDeriv::grad - parameters must be the first argument
  # and it must return a scalar or vector (hrf_lwu returns a vector of length(t_val))
  target_func_for_grad <- function(params_vec, t_val) {
    hrf_lwu(t = t_val, tau = params_vec[1], sigma = params_vec[2], rho = params_vec[3], normalize = "none")
  }

  # Calculate h0 (the HRF at theta0)
  h0 <- hrf_lwu(t = t, tau = theta0["tau"], sigma = theta0["sigma"], rho = theta0["rho"], normalize = "none")

  if (normalize_primary == "height") {
    max_abs_h0 <- max(abs(h0), na.rm = TRUE)
    if (max_abs_h0 > 1e-10) {
      h0 <- h0 / max_abs_h0
    }
  }

  # Calculate partial derivatives using numDeriv::grad
  # grad() will iterate over each element of `t` if `target_func_for_grad` is vectorized over t,
  # which it is. We want the gradient for each time point.
  # However, numDeriv::grad expects the function to return a single scalar for jacobian calculation.
  # So, we must loop over t points for numDeriv.

  deriv_matrix <- matrix(NA, nrow = length(t), ncol = 3)
  colnames(deriv_matrix) <- c("d_tau", "d_sigma", "d_rho")

  for (i in seq_along(t)) {
    # numDeriv::grad needs a function that takes params and returns a SINGLE value
    # So we create a wrapper for each time point t[i]
    current_t_func <- function(params_vec) {
      hrf_lwu(t = t[i], tau = params_vec[1], sigma = params_vec[2], rho = params_vec[3], normalize = "none")
    }
    # Calculate gradient (vector of 3 partial derivatives) at t[i] w.r.t theta0
    grad_at_t_i <- numDeriv::grad(func = current_t_func, x = theta0)
    deriv_matrix[i, ] <- grad_at_t_i
  }

  basis_mat <- cbind(h0 = h0, deriv_matrix)
  return(basis_mat)
}
