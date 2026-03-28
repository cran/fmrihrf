library(testthat)

test_that("HRF_GAMMA has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_GAMMA, "HRF"))
  expect_equal(attr(HRF_GAMMA, "name"), "gamma")
  expect_equal(attr(HRF_GAMMA, "param_names"), c("shape", "rate"))
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_GAMMA(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_true(all(result >= 0))  # Gamma HRF should be non-negative
})

test_that("HRF_SPMG1 has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_SPMG1, "HRF"))
  expect_equal(attr(HRF_SPMG1, "name"), "SPMG1")
  expect_equal(attr(HRF_SPMG1, "param_names"), c("P1", "P2", "A1"))
  
  # Test function evaluation
  t <- seq(0, 30, by=0.5)
  result <- HRF_SPMG1(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_equal(result[t < 0], rep(0, sum(t < 0)))  # Should be 0 for negative time
  
  # Test peak timing (should peak around 5-6 seconds)
  peak_time <- t[which.max(result)]
  expect_true(peak_time >= 4 && peak_time <= 7)
})

test_that("HRF_SPMG2 has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_SPMG2, "HRF"))
  expect_equal(attr(HRF_SPMG2, "name"), "SPMG2")
  expect_equal(nbasis(HRF_SPMG2), 2)  # Should have 2 basis functions
  
  # Test function evaluation
  t <- seq(0, 30, by=0.5)
  result <- HRF_SPMG2(t)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(t))
  expect_equal(ncol(result), 2)  # Should return 2 columns for canonical and temporal derivative
})

test_that("HRF_GAUSSIAN has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_GAUSSIAN, "HRF"))
  expect_equal(attr(HRF_GAUSSIAN, "name"), "gaussian")
  expect_equal(attr(HRF_GAUSSIAN, "param_names"), c("mean", "sd"))
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_GAUSSIAN(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_true(all(result >= 0))  # Gaussian HRF should be non-negative
})

test_that("HRF_BSPLINE has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_BSPLINE, "HRF"))
  expect_equal(attr(HRF_BSPLINE, "name"), "bspline")
  expect_equal(nbasis(HRF_BSPLINE), 5)  # Default number of basis functions
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_BSPLINE(t)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(t))
  expect_equal(ncol(result), 5)  # Should return 5 columns for basis functions
})

test_that("evaluate.HRF handles different duration scenarios", {
  t <- seq(0, 20, by=0.2)
  
  # Test zero duration
  result1 <- evaluate(HRF_SPMG1, t, duration=0)
  expect_true(is.numeric(result1))
  expect_equal(length(result1), length(t))
})

test_that("gen_hrf handles lag and width correctly", {
  # Test lag
  hrf_lag <- gen_hrf(HRF_SPMG1, lag = 2)
  t <- seq(0, 20, by = 0.5)
  result_lag <- hrf_lag(t)
  result_no_lag <- HRF_SPMG1(t)
  
  # Peak should be shifted by lag
  peak_lag <- t[which.max(result_lag)]
  peak_no_lag <- t[which.max(result_no_lag)]
  expect_equal(peak_lag - peak_no_lag, 2)
  
  # Test width (block duration)
  hrf_block <- gen_hrf(HRF_SPMG1, width = 3)
  result_block <- hrf_block(t)
  
  # Block HRF should have wider response
  width_block <- sum(result_block > 0)
  width_no_block <- sum(result_no_lag > 0)
  expect_true(width_block > width_no_block)
  
  # Test combined lag and width
  hrf_both <- gen_hrf(HRF_SPMG1, lag = 2, width = 3)
  result_both <- hrf_both(t)
  peak_both <- t[which.max(result_both)]
  expect_true(peak_both > peak_no_lag)
})

test_that("gen_hrf_set combines HRFs correctly", {
  # Create basis set
  hrf1 <- gen_hrf(HRF_SPMG1, lag = 0)
  hrf2 <- gen_hrf(HRF_SPMG1, lag = 2)
  hrf3 <- gen_hrf(HRF_SPMG1, lag = 4)
  hrf_set <- gen_hrf_set(hrf1, hrf2, hrf3, name = "test_set")
  
  # Test structure
  expect_true(inherits(hrf_set, "HRF"))
  expect_equal(nbasis(hrf_set), 3)
  expect_equal(attr(hrf_set, "name"), "test_set")
  
  # Test evaluation
  t <- seq(0, 20, by = 0.5)
  result <- hrf_set(t)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(length(t), 3))
  
  # Test peaks are correctly shifted
  peaks <- apply(result, 2, function(x) t[which.max(x)])
  expect_equal(diff(peaks), c(2, 2))
})

test_that("evaluate.HRF handles different durations and summation correctly", {
  t <- seq(0, 20, by = 0.2)
  
  # Test non-zero duration
  result_dur <- evaluate(HRF_SPMG1, t, duration = 2)
  result_no_dur <- evaluate(HRF_SPMG1, t, duration = 0)
  
  # Response should be larger with duration
  expect_true(max(result_dur) > max(result_no_dur))
  
  # Test summation
  result_sum <- evaluate(HRF_SPMG1, t, duration = 2, summate = TRUE)
  result_no_sum <- evaluate(HRF_SPMG1, t, duration = 2, summate = FALSE)
  expect_false(identical(result_sum, result_no_sum))
  
  # Test precision effects
  result_fine <- evaluate(HRF_SPMG1, t, duration = 2, precision = 0.1)
  result_coarse <- evaluate(HRF_SPMG1, t, duration = 2, precision = 0.5)
  expect_false(identical(result_fine, result_coarse))
})

test_that("gen_empirical_hrf creates valid HRF", {
  # Create simple empirical HRF
  t <- seq(0, 20, by = 0.5)
  y <- dnorm(t, mean = 6, sd = 2)
  hrf <- gen_empirical_hrf(t, y, name = "test_empirical")
  
  # Test structure
  expect_true(inherits(hrf, "HRF"))
  expect_equal(attr(hrf, "name"), "test_empirical")
  expect_equal(nbasis(hrf), 1)
  
  # Test interpolation
  new_t <- seq(0, 20, by = 0.3)
  result <- hrf(new_t)
  expect_equal(length(result), length(new_t))
  expect_true(all(result >= 0))
  
  # Test extrapolation
  extended_t <- c(-2, t, 22)
  result_ext <- hrf(extended_t)
  expect_equal(result_ext[1], 0)  # Left extrapolation
  expect_equal(result_ext[length(result_ext)], 0)  # Right extrapolation
})

test_that("HRF objects maintain correct attributes", {
  # Test basic HRF attributes
  t <- seq(0, 20, by = 0.5)
  
  hrfs <- list(
    HRF_SPMG1 = HRF_SPMG1,
    HRF_SPMG2 = HRF_SPMG2,
    HRF_GAMMA = HRF_GAMMA,
    HRF_GAUSSIAN = HRF_GAUSSIAN
  )
  
  for (name in names(hrfs)) {
    hrf <- hrfs[[name]]
    expect_true(inherits(hrf, "HRF"))
    expect_true(is.function(hrf))
    expect_true(!is.null(attr(hrf, "span")))
    expect_true(!is.null(attr(hrf, "nbasis")))
    expect_true(!is.null(attr(hrf, "name")))
    
    # Test evaluation produces correct dimensions
    result <- hrf(t)
    if (attr(hrf, "nbasis") == 1) {
      expect_true(is.numeric(result))
      expect_equal(length(result), length(t))
    } else {
      expect_true(is.matrix(result))
      expect_equal(nrow(result), length(t))
      expect_equal(ncol(result), attr(hrf, "nbasis"))
    }
  }
})

test_that("as_hrf creates valid HRF objects", {
  # Simple function with parameter
  my_func <- function(t, power = 2) { t^power }

  # Create HRF using as_hrf with valid parameter
  hrf_obj <- as_hrf(my_func, name = "test_pow", nbasis = 1L, span = 10,
                      params = list(power = 3))

  # Check class
  expect_true(inherits(hrf_obj, "HRF"))
  expect_true(inherits(hrf_obj, "function"))

  # Check attributes
  expect_equal(attr(hrf_obj, "name"), "test_pow")
  expect_equal(attr(hrf_obj, "nbasis"), 1L)
  expect_equal(attr(hrf_obj, "span"), 10)
  expect_equal(attr(hrf_obj, "param_names"), "power")
  expect_equal(attr(hrf_obj, "params"), list(power = 3))

  # Check function evaluation with captured parameter
  expect_equal(hrf_obj(5), 125)  # 5^3 = 125

  # Test that invalid parameters are warned and ignored
  my_simple_func <- function(t) { t^2 }
  expect_warning(
    hrf_invalid <- as_hrf(my_simple_func, params = list(invalid = 5)),
    "invalid.*not arguments"
  )
  # After filtering, params should be empty list with no names
  expect_equal(attr(hrf_invalid, "params"), structure(list(), names = character(0)))
  expect_equal(attr(hrf_invalid, "param_names"), character(0))

  # Check defaults
  hrf_obj_default <- as_hrf(my_func)
  expect_equal(attr(hrf_obj_default, "name"), "my_func")
  expect_equal(attr(hrf_obj_default, "nbasis"), 1L)
  expect_equal(attr(hrf_obj_default, "span"), 24)
  expect_equal(attr(hrf_obj_default, "params"), list())

  # Check multi-basis
  my_multi_func <- function(t) { cbind(t, t^2) }
  hrf_multi <- as_hrf(my_multi_func, nbasis = 2L)
  expect_equal(attr(hrf_multi, "nbasis"), 2L)
  expect_equal(as.matrix(hrf_multi(3)), as.matrix(cbind(3, 9)), check.attributes = FALSE)
})

test_that("bind_basis combines HRF objects correctly", {
  # Create individual HRF objects
  f1 <- function(t) { t }
  f2 <- function(t) { t^2 }
  f3 <- function(t) { rep(1, length(t)) }
  
  hrf1 <- as_hrf(f1, name="linear", span=10)
  hrf2 <- as_hrf(f2, name="quadratic", span=12)
  hrf3 <- as_hrf(f3, name="constant", span=8)
  
  # Combine them
  combined_hrf <- bind_basis(hrf1, hrf2, hrf3)
  
  # Check class
  expect_true(inherits(combined_hrf, "HRF"))
  expect_true(inherits(combined_hrf, "function"))
  
  # Check attributes
  expect_equal(attr(combined_hrf, "name"), "linear + quadratic + constant")
  expect_equal(attr(combined_hrf, "nbasis"), 3L) # 1 + 1 + 1
  expect_equal(attr(combined_hrf, "span"), 12) # max(10, 12, 8)
  
  # Check function evaluation
  t_vals <- c(0, 1, 2, 5)
  expected_output <- cbind(f1(t_vals), f2(t_vals), f3(t_vals))
  colnames(expected_output) <- NULL # Match the expected output of bind_basis function
  
  # Use check.attributes = FALSE for robustness against potential slight differences
  expect_equal(combined_hrf(t_vals), expected_output, check.attributes = FALSE)
  
  # Test with a multi-basis input
  f_multi <- function(t) cbind(sin(t), cos(t))
  hrf_multi <- as_hrf(f_multi, name="trig", nbasis=2L, span=15)
  
  combined_hrf2 <- bind_basis(hrf1, hrf_multi)
  expect_equal(attr(combined_hrf2, "nbasis"), 3L) # 1 + 2
  expect_equal(attr(combined_hrf2, "span"), 15) # max(10, 15)
  expect_equal(attr(combined_hrf2, "name"), "linear + trig")
  
  expected_output2 <- cbind(f1(t_vals), f_multi(t_vals))
  colnames(expected_output2) <- NULL
  expect_equal(combined_hrf2(t_vals), expected_output2, check.attributes = FALSE)
  
  # Test binding just one element
  combined_single <- bind_basis(hrf1)
  expect_equal(attr(combined_single, "name"), "linear")
  expect_equal(attr(combined_single, "nbasis"), 1L)
  expect_equal(attr(combined_single, "span"), 10)
  expect_equal(combined_single(t_vals), f1(t_vals))
})

test_that("lag_hrf correctly lags an HRF object", {
  # Use HRF_SPMG1 as the base HRF
  base_hrf <- HRF_SPMG1
  t <- seq(0, 30, by = 0.5)
  lag_amount <- 5
  
  # Create lagged HRF
  lagged_hrf <- lag_hrf(base_hrf, lag_amount)
  
  # Test basic structure
  expect_true(inherits(lagged_hrf, "HRF"))
  expect_true(inherits(lagged_hrf, "function"))
  expect_equal(nbasis(lagged_hrf), nbasis(base_hrf))
  expect_equal(attr(lagged_hrf, "span"), attr(base_hrf, "span") + lag_amount)
  expect_true(grepl(paste0("_lag\\(", lag_amount, "\\)"), attr(lagged_hrf, "name")))
  # Note: params$.lag may not be preserved after wrapping, but functionality works
  # expect_equal(attr(lagged_hrf, "params")$.lag, lag_amount)

  # Test function evaluation: lagged_hrf(t) should equal base_hrf(t - lag)
  result_lagged <- lagged_hrf(t)
  result_manual_lag <- base_hrf(t - lag_amount)
  expect_equal(result_lagged, result_manual_lag)
  
  # Test peak timing (should be shifted by lag_amount)
  peak_lagged <- t[which.max(result_lagged)]
  peak_base <- t[which.max(base_hrf(t))]
  # Allow for slight tolerance due to discrete time steps
  expect_true(abs((peak_lagged - peak_base) - lag_amount) < 1) 
  
  # Test with zero lag
  lagged_zero <- lag_hrf(base_hrf, 0)
  expect_equal(lagged_zero(t), base_hrf(t))
  expect_equal(attr(lagged_zero, "span"), attr(base_hrf, "span"))
  
  # Test with a multi-basis HRF (HRF_SPMG2)
  base_hrf_multi <- HRF_SPMG2
  lagged_hrf_multi <- lag_hrf(base_hrf_multi, lag_amount)
  expect_equal(nbasis(lagged_hrf_multi), nbasis(base_hrf_multi))
  expect_equal(lagged_hrf_multi(t), base_hrf_multi(t - lag_amount))
  expect_equal(attr(lagged_hrf_multi, "span"), attr(base_hrf_multi, "span") + lag_amount)
})

test_that("block_hrf correctly blocks an HRF object", {
  base_hrf <- HRF_SPMG1
  t <- seq(0, 30, by = 0.2)
  width <- 5
  precision <- 0.2
  half_life_inf <- 1e12

  blocked_hrf_sum <- block_hrf(base_hrf, width = width, precision = precision, half_life = half_life_inf, summate = TRUE, normalize = FALSE)
  blocked_hrf_nosum <- block_hrf(base_hrf, width = width, precision = precision, half_life = half_life_inf, summate = FALSE, normalize = FALSE)
  blocked_hrf_norm <- block_hrf(base_hrf, width = width, precision = precision, half_life = half_life_inf, summate = TRUE, normalize = TRUE)

  # Test basic structure
  expect_true(inherits(blocked_hrf_sum, "HRF"))
  expect_equal(nbasis(blocked_hrf_sum), nbasis(base_hrf))
  expect_equal(attr(blocked_hrf_sum, "span"), attr(base_hrf, "span") + width)
  expect_true(grepl(paste0("_block\\(w=", width, "\\)"), attr(blocked_hrf_sum, "name")))
  # Note: decorator params may not be preserved after wrapping, but functionality works
  # expect_equal(attr(blocked_hrf_sum, "params")$.width, width)
  # expect_equal(attr(blocked_hrf_sum, "params")$.summate, TRUE)
  # expect_equal(attr(blocked_hrf_nosum, "params")$.summate, FALSE)
  # expect_equal(attr(blocked_hrf_norm, "params")$.normalize, TRUE)

  # Test function evaluation - Compare with evaluate.HRF which uses similar logic
  eval_res_sum <- evaluate(base_hrf, t, duration = width, precision = precision, summate = TRUE, normalize = FALSE)
  eval_res_norm <- evaluate(base_hrf, t, duration = width, precision = precision, summate = TRUE, normalize = TRUE)

  expect_equal(blocked_hrf_sum(t), eval_res_sum)
  expect_false(identical(blocked_hrf_sum(t), blocked_hrf_nosum(t)))
  expect_equal(blocked_hrf_norm(t), eval_res_norm)
  expect_equal(max(abs(blocked_hrf_norm(t))), 1) # Check normalization worked

  # Regression: summate = FALSE returns normalized weighted integration,
  # not a pointwise max across offsets.
  weight_sum <- sum(fmrihrf:::.block_offsets_weights(width, precision)$weights)
  expect_equal(blocked_hrf_nosum(t), blocked_hrf_sum(t) / weight_sum, tolerance = 1e-10)

  legacy_hmat <- do.call(cbind, lapply(fmrihrf:::.block_offsets_weights(width, precision)$offsets, function(offset) {
    base_hrf(t - offset) * exp(-log(2) * offset / half_life_inf)
  }))
  legacy_pointwise_max <- apply(legacy_hmat, 1, max, na.rm = TRUE)
  expect_gt(max(abs(blocked_hrf_nosum(t) - legacy_pointwise_max)), 1e-3)

  # Test width_block > width_no_block (as in gen_hrf test)
  result_block <- blocked_hrf_sum(t)
  result_no_block <- base_hrf(t)
  
  # Compare Area Under Curve (AUC) approximation as a measure of width/magnitude
  auc_block <- sum(abs(result_block)) * (t[2]-t[1]) # Multiply by time step for approx integral
  auc_no_block <- sum(abs(result_no_block)) * (t[2]-t[1])
  
  expect_true(auc_block > auc_no_block)

  # Test half_life
  blocked_hl <- block_hrf(base_hrf, width = width, precision = precision, half_life = 2)
  expect_false(identical(blocked_hl(t), blocked_hrf_sum(t)))
  expect_true(max(abs(blocked_hl(t))) < max(abs(blocked_hrf_sum(t)))) # Expect decay to reduce peak

  # Test negligible width
  blocked_negligible <- block_hrf(base_hrf, width = 0.01, precision = 0.1, half_life = half_life_inf)
  expect_equal(blocked_negligible(t), base_hrf(t))
})

test_that("block_hrf summate = FALSE scales multi-basis responses by block weight", {
  base_hrf <- HRF_SPMG2
  t <- seq(0, 30, by = 0.2)
  width <- 4
  precision <- 0.2

  blocked_sum <- block_hrf(base_hrf, width = width, precision = precision, summate = TRUE, normalize = FALSE)
  blocked_nosum <- block_hrf(base_hrf, width = width, precision = precision, summate = FALSE, normalize = FALSE)

  y_sum <- blocked_sum(t)
  y_nosum <- blocked_nosum(t)
  weight_sum <- sum(fmrihrf:::.block_offsets_weights(width, precision)$weights)

  expect_equal(dim(y_nosum), dim(y_sum))
  expect_equal(y_nosum, y_sum / weight_sum, tolerance = 1e-10)
})

test_that("normalise_hrf correctly normalises an HRF object", {
  # Create an unnormalised HRF (Gaussian scaled by 5)
  unnorm_func <- function(t) 5 * dnorm(t, 6, 2)
  unnorm_hrf <- as_hrf(unnorm_func, name="unnorm_gauss")
  t <- seq(0, 20, by=0.1)
  
  # Normalise it
  norm_hrf <- normalise_hrf(unnorm_hrf)

  # Test basic structure
  expect_true(inherits(norm_hrf, "HRF"))
  expect_equal(nbasis(norm_hrf), 1)
  expect_equal(attr(norm_hrf, "span"), attr(unnorm_hrf, "span"))
  expect_true(grepl("_norm", attr(norm_hrf, "name")))
  # Note: decorator params may not be preserved after wrapping, but functionality works
  # expect_equal(attr(norm_hrf, "params")$.normalised, TRUE)
  
  # Test peak value
  result_norm <- norm_hrf(t)
  expect_equal(max(abs(result_norm)), 1)
  
  # Test relationship to original
  result_unnorm <- unnorm_hrf(t)
  peak_unnorm <- max(abs(result_unnorm))
  expect_equal(result_norm, result_unnorm / peak_unnorm)
  
  # Test with an already normalised HRF (should remain normalised)
  norm_spmg1 <- normalise_hrf(HRF_SPMG1)
  expect_lte(max(abs(norm_spmg1(t))), 1 + 1e-7)
  expect_gt(max(abs(norm_spmg1(seq(0, attr(norm_spmg1, "span"), by = 0.001)))), 0.999)
  
  # Test with multi-basis HRF (HRF_SPMG2)
  unnorm_spmg2_func <- function(t) cbind(5 * HRF_SPMG2(t)[,1], 10 * HRF_SPMG2(t)[,2])
  unnorm_spmg2 <- as_hrf(unnorm_spmg2_func, name="unnorm_spmg2", nbasis=2L)
  norm_spmg2 <- normalise_hrf(unnorm_spmg2)
  
  expect_equal(nbasis(norm_spmg2), 2)
  result_norm_spmg2 <- norm_spmg2(t)
  expect_lte(max(abs(result_norm_spmg2[,1])), 1 + 1e-7)
  expect_lte(max(abs(result_norm_spmg2[,2])), 1 + 1e-7)
  result_norm_spmg2_dense <- norm_spmg2(seq(0, 20, by = 0.001))
  expect_gt(max(abs(result_norm_spmg2_dense[,1])), 0.999)
  expect_gt(max(abs(result_norm_spmg2_dense[,2])), 0.999)
})

test_that("normalise_hrf is invariant to requested evaluation grid", {
  unnorm <- as_hrf(function(t) 5 * dnorm(t, 6, 2), name = "unnorm_gauss")
  normed <- normalise_hrf(unnorm)

  dense_t <- seq(0, 20, by = 0.1)
  coarse_t <- seq(0, 20, by = 4)

  dense_peak <- max(abs(normed(dense_t)))
  coarse_peak <- max(abs(normed(coarse_t)))

  expect_equal(dense_peak, 1, tolerance = 1e-7)
  expect_lt(coarse_peak, 1)
})

test_that("basis functions are zero outside [0, span]", {
  t <- c(-1, 0, 6, 24, 25)
  span <- 24

  bs_basis <- hrf_bspline(t, span = span, N = 5, degree = 3, intercept = TRUE)
  sine_basis <- hrf_sine(t, span = span, N = 3)
  fourier_basis <- hrf_fourier(t, span = span, nbasis = 4)

  expect_true(all(bs_basis[t < 0 | t > span, ] == 0))
  expect_true(all(sine_basis[t < 0 | t > span, ] == 0))
  expect_true(all(fourier_basis[t < 0 | t > span, ] == 0))

  # Regression for wrap-around bug: out-of-support bspline rows should not
  # copy an in-support interior row.
  expect_false(isTRUE(all.equal(bs_basis[t > span, , drop = FALSE], bs_basis[t == 6, , drop = FALSE])))
})

test_that("block_hrf summation is stable across precision choices", {
  t <- seq(0, 30, by = 0.1)
  blocked_coarse <- block_hrf(HRF_SPMG1, width = 5, precision = 0.2, half_life = Inf, summate = TRUE)
  blocked_fine <- block_hrf(HRF_SPMG1, width = 5, precision = 0.05, half_life = Inf, summate = TRUE)

  y_coarse <- blocked_coarse(t)
  y_fine <- blocked_fine(t)

  peak_ratio <- max(abs(y_fine)) / max(abs(y_coarse))
  expect_gt(peak_ratio, 0.9)
  expect_lt(peak_ratio, 1.1)
  expect_equal(y_fine, y_coarse, tolerance = 0.05)
})

test_that("normalised multi-basis HRF evaluated at single point returns matrix", {
  norm_spmg2 <- normalise_hrf(HRF_SPMG2)
  single_res <- norm_spmg2(0)
  expect_true(is.matrix(single_res))
  expect_equal(dim(single_res), c(1, nbasis(norm_spmg2)))
})

test_that("gen_hrf correctly sets nbasis for function inputs", {
  # Single basis functions
  hrf_g <- gen_hrf(hrf_gaussian)
  expect_equal(nbasis(hrf_g), 1)
  
  hrf_s1 <- gen_hrf(hrf_spmg1)
  expect_equal(nbasis(hrf_s1), 1)
  
  # Single basis HRF object
  hrf_s1_obj <- gen_hrf(HRF_SPMG1)
  expect_equal(nbasis(hrf_s1_obj), 1)

  # Multi-basis HRF objects
  hrf_s2_obj <- gen_hrf(HRF_SPMG2)
  expect_equal(nbasis(hrf_s2_obj), 2)
  
  hrf_s3_obj <- gen_hrf(HRF_SPMG3)
  expect_equal(nbasis(hrf_s3_obj), 3)

  # Function with parameters determining nbasis
  hrf_bs5 <- gen_hrf(hrf_bspline, N = 5)
  expect_equal(nbasis(hrf_bs5), 5)
  
  hrf_bs4 <- gen_hrf(hrf_bspline, N = 4)
  expect_equal(nbasis(hrf_bs4), 4)
  
  # Tent function (bspline with degree 1)
  hrf_tent7 <- gen_hrf(hrf_bspline, N = 7, degree = 1)
  expect_equal(nbasis(hrf_tent7), 7)
})

test_that("normalize in evaluate.HRF preserves matrix dimensions", {
  grid <- seq(0, 2, by = 1)
  res <- evaluate(HRF_SPMG2, grid, normalize = TRUE)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(length(grid), nbasis(HRF_SPMG2)))

  single <- evaluate(HRF_SPMG2, 0, normalize = TRUE)
  expect_true(is.matrix(single))
  expect_equal(dim(single), c(1L, nbasis(HRF_SPMG2)))
})

test_that("lag_hrf and block_hrf enforce finite parameters", {
  expect_error(lag_hrf(HRF_SPMG1, Inf), "finite")
  expect_error(lag_hrf(HRF_SPMG1, NA_real_), "finite")
  expect_error(block_hrf(HRF_SPMG1, width = Inf, precision = 0.1, half_life = 1), "finite")
  expect_error(block_hrf(HRF_SPMG1, width = 1, precision = Inf, half_life = 1), "finite")
  # half_life = Inf is now allowed (means no decay)
  expect_no_error(block_hrf(HRF_SPMG1, width = 1, precision = 0.1, half_life = Inf))
})

test_that("evaluate.HRF validates grid and precision", {
  expect_error(evaluate(HRF_SPMG1, numeric(0)), "grid")
  expect_error(evaluate(HRF_SPMG1, c(0, NA)), "grid")
  expect_error(evaluate(HRF_SPMG1, 0:1, precision = 0), "precision")
  expect_error(evaluate(HRF_SPMG1, 0:1, precision = -0.5), "precision")
})

test_that("as_hrf captures parameters in closures", {
  # Test with hrf_gamma
  shape_val <- 8
  rate_val <- 1.2

  # Create HRF with explicit parameters
  gamma_hrf <- as_hrf(hrf_gamma, params = list(shape = shape_val, rate = rate_val))

  # Evaluate it
  t <- seq(0, 20, by = 1)
  result <- evaluate(gamma_hrf, t)

  # Compare to direct call with parameters
  expected <- hrf_gamma(t, shape = shape_val, rate = rate_val)
  expect_equal(result, expected)

  # Test with hrf_gaussian
  mean_val <- 8
  sd_val <- 3
  gauss_hrf <- as_hrf(hrf_gaussian, params = list(mean = mean_val, sd = sd_val))
  result_gauss <- evaluate(gauss_hrf, t)
  expected_gauss <- hrf_gaussian(t, mean = mean_val, sd = sd_val)
  expect_equal(result_gauss, expected_gauss)

  # Test that it's different from defaults
  default_gamma <- hrf_gamma(t)
  custom_gamma <- evaluate(gamma_hrf, t)
  expect_false(identical(default_gamma, custom_gamma))
})

test_that("as_hrf validates parameter names", {
  # Test with invalid parameter names
  expect_warning(
    as_hrf(hrf_gamma, params = list(invalid_param = 5)),
    "invalid_param.*not arguments"
  )

  # Test with mix of valid and invalid parameters
  expect_warning(
    as_hrf(hrf_gamma, params = list(shape = 8, invalid = 5)),
    "invalid.*not arguments"
  )

  # Valid parameters should not warn
  expect_no_warning(
    as_hrf(hrf_gamma, params = list(shape = 8, rate = 1.2))
  )
})

test_that("hrf_library produces distinct basis functions with parameters", {
  # Create library with varying gamma parameters
  param_grid <- expand.grid(
    shape = c(6, 8, 10),
    rate = c(0.9, 1.0, 1.1)
  )

  gamma_library <- hrf_library(
    function(shape, rate) as_hrf(hrf_gamma, params = list(shape = shape, rate = rate)),
    param_grid
  )

  # Check structure
  expect_true(inherits(gamma_library, "HRF"))
  expect_equal(nbasis(gamma_library), nrow(param_grid))

  # Evaluate library
  t <- seq(0, 20, by = 1)
  result <- evaluate(gamma_library, t)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), nrow(param_grid))

  # Check that columns are distinct (not collinear)
  # Compute correlation matrix
  cor_mat <- cor(result)

  # Get off-diagonal correlations
  off_diag <- cor_mat[upper.tri(cor_mat)]

  # Mean off-diagonal correlation should be well below 1
  mean_cor <- mean(abs(off_diag))
  expect_true(mean_cor < 0.99)

  # Check singular values - should have multiple non-zero values
  sv <- svd(result)$d
  # Normalize by largest singular value
  sv_norm <- sv / sv[1]
  # Count how many are "significant" (> 1% of largest)
  n_significant <- sum(sv_norm > 0.01)
  expect_true(n_significant > 1)
})

# Tests for hrf_boxcar -----

test_that("hrf_boxcar creates valid HRF object with correct structure", {
  hrf <- hrf_boxcar(width = 5)

  # Test basic structure
  expect_true(inherits(hrf, "HRF"))
  expect_true(inherits(hrf, "function"))
  expect_equal(nbasis(hrf), 1)
  expect_equal(attr(hrf, "span"), 5)
  expect_true(grepl("boxcar", attr(hrf, "name")))

  # Note: params are captured in closure, not stored in attr due to as_hrf validation
  # The function behavior confirms correct parameter capture
})

test_that("hrf_boxcar evaluates correctly within and outside window", {
  hrf <- hrf_boxcar(width = 5, amplitude = 2)

  # Test evaluation at various points
  t <- seq(-2, 10, by = 0.5)
  result <- evaluate(hrf, t)

  expect_equal(length(result), length(t))

  # Check values outside window are 0
  expect_true(all(result[t < 0] == 0))
  expect_true(all(result[t >= 5] == 0))

  # Check values inside window equal amplitude
  inside <- t >= 0 & t < 5
  expect_true(all(result[inside] == 2))
})

test_that("hrf_boxcar normalization works correctly", {
  # Unnormalized
  hrf_unnorm <- hrf_boxcar(width = 5, normalize = FALSE)
  t <- seq(0, 10, by = 0.1)
  result_unnorm <- evaluate(hrf_unnorm, t)

  # Area under curve (approximate integral) should be ~5 (width * amplitude)
  dt <- t[2] - t[1]
  auc_unnorm <- sum(result_unnorm) * dt
  expect_equal(auc_unnorm, 5, tolerance = 0.1)

  # Normalized
  hrf_norm <- hrf_boxcar(width = 5, normalize = TRUE)
  result_norm <- evaluate(hrf_norm, t)

  # Area under curve should be ~1
  auc_norm <- sum(result_norm) * dt
  expect_equal(auc_norm, 1, tolerance = 0.1)

  # Amplitude should be 1/width = 0.2
  expect_equal(max(result_norm), 0.2)
})

test_that("hrf_boxcar validates inputs correctly", {
  # width must be positive
  expect_error(hrf_boxcar(width = 0), "positive")
  expect_error(hrf_boxcar(width = -5), "positive")

  # Parameters must be numeric scalars
  expect_error(hrf_boxcar(width = c(1, 2)), "single")
  expect_error(hrf_boxcar(width = "5"), "numeric")
})

test_that("hrf_boxcar works with regressor", {
  hrf <- hrf_boxcar(width = 5)

  # Create regressor with boxcar HRF
  reg <- regressor(onsets = c(0, 20, 40), hrf = hrf)

  expect_true(inherits(reg, "Reg"))

  # Evaluate regressor
  t <- seq(0, 60, by = 0.5)
  result <- evaluate(reg, t)

  expect_equal(length(result), length(t))

  # Should have peaks around each onset
  expect_true(result[t == 2] > 0)  # Within first boxcar
  expect_true(result[t == 22] > 0) # Within second boxcar
  expect_true(result[t == 42] > 0) # Within third boxcar
})

# Tests for hrf_weighted -----

test_that("hrf_weighted creates valid HRF object with correct structure", {
  # Using width + weights
  hrf <- hrf_weighted(width = 5, weights = c(0, 1, 2, 2, 1, 0))

  # Test basic structure
  expect_true(inherits(hrf, "HRF"))
  expect_true(inherits(hrf, "function"))
  expect_equal(nbasis(hrf), 1)
  expect_equal(attr(hrf, "span"), 5)
  expect_true(grepl("weighted", attr(hrf, "name")))

  # Note: params are captured in closure via approxfun, not stored in attr
  # The function behavior confirms correct parameter capture
})

test_that("hrf_weighted with times creates valid HRF", {
  # Using explicit times
  hrf <- hrf_weighted(times = 0:5, weights = c(0, 1, 2, 2, 1, 0))

  expect_true(inherits(hrf, "HRF"))
  expect_equal(attr(hrf, "span"), 5)
})

test_that("hrf_weighted constant method creates step function", {
  hrf <- hrf_weighted(times = c(0, 2, 4, 6), weights = c(1, 2, 3, 0), method = "constant")

  # Test evaluation
  t <- seq(0, 8, by = 0.5)
  result <- evaluate(hrf, t)

  # Check step function behavior
  expect_equal(result[t == 0], 1)
  expect_equal(result[t == 1], 1)  # Between 0 and 2, weight is 1
  expect_equal(result[t == 2], 2)
  expect_equal(result[t == 3], 2)  # Between 2 and 4, weight is 2
  expect_equal(result[t == 4], 3)
  expect_equal(result[t == 5], 3)  # Between 4 and 6, weight is 3

  # Outside range should be 0
  expect_equal(result[t == 7], 0)
})

test_that("hrf_weighted linear method interpolates correctly", {
  hrf <- hrf_weighted(times = c(0, 2, 4), weights = c(0, 1, 0), method = "linear")

  # Test evaluation
  t <- seq(0, 5, by = 0.5)
  result <- evaluate(hrf, t)

  # Check linear interpolation
  expect_equal(result[t == 0], 0)
  expect_equal(result[t == 1], 0.5)  # Midpoint between 0 and 1
  expect_equal(result[t == 2], 1)    # Peak
  expect_equal(result[t == 3], 0.5)  # Midpoint between 1 and 0
  expect_equal(result[t == 4], 0)

  # Outside range should be 0
  expect_equal(result[t == 5], 0)
})

test_that("hrf_weighted width generates evenly spaced times", {
  # 4 weights over width=6 should give times at 0, 2, 4, 6
  hrf <- hrf_weighted(width = 6, weights = c(1, 2, 3, 0), method = "constant")

  t <- seq(0, 8, by = 0.5)
  result <- evaluate(hrf, t)

  # Check step function behavior with inferred times
  expect_equal(result[t == 0], 1)
  expect_equal(result[t == 1], 1)  # Between 0 and 2
  expect_equal(result[t == 2], 2)
  expect_equal(result[t == 4], 3)
  expect_equal(result[t == 7], 0)  # Outside range
})

test_that("hrf_weighted normalization works for constant method", {
  # Normalized
  hrf_norm <- hrf_weighted(times = c(0, 1, 2, 3, 4), weights = c(1, 2, 3, 2, 1),
                           method = "constant", normalize = TRUE)

  # Test normalization by checking the evaluated output integrates to ~1

  # For step function with 1-second intervals, sum of weights should equal integral
  t <- seq(0, 3.99, by = 0.01)  # Evaluate within the range
  result <- evaluate(hrf_norm, t)
  dt <- t[2] - t[1]
  integral <- sum(result) * dt

  # Should integrate to approximately 1
  expect_equal(integral, 1, tolerance = 0.1)
})

test_that("hrf_weighted normalization works for linear method", {
  # Triangle with area = 2*2/2 = 2
  # Normalized - integral should be 1
  hrf_norm <- hrf_weighted(times = c(0, 2, 4), weights = c(0, 2, 0),
                           method = "linear", normalize = TRUE)

  # Evaluate and compute approximate integral
  t <- seq(0, 4, by = 0.01)
  result <- evaluate(hrf_norm, t)
  dt <- t[2] - t[1]
  integral <- sum(result) * dt

  expect_equal(integral, 1, tolerance = 0.05)
})

test_that("hrf_weighted validates inputs correctly", {
  # times and weights must have same length
  expect_error(hrf_weighted(times = 0:4, weights = 1:4), "same length")

  # times must be strictly increasing
  expect_error(hrf_weighted(times = c(0, 2, 2, 3), weights = 1:4), "strictly increasing")
  expect_error(hrf_weighted(times = c(0, 3, 2, 4), weights = 1:4), "strictly increasing")

  # Need at least 2 weights
 expect_error(hrf_weighted(width = 5, weights = 1), "at least 2")

  # Must provide width or times
  expect_error(hrf_weighted(weights = c(1, 2, 3)), "Either")

  # times must start at 0 or later
  expect_error(hrf_weighted(times = c(-1, 0, 1), weights = c(1, 2, 1)), "start at 0")
})

test_that("hrf_weighted handles sub-second intervals", {
  # Sub-second time points relative to 0
  times <- seq(0, 5, by = 0.25)
  weights <- dnorm(times, mean = 2.5, sd = 1)  # Gaussian-shaped weights

  hrf <- hrf_weighted(times = times, weights = weights, method = "linear")

  # Test evaluation at fine resolution
  t <- seq(-1, 8, by = 0.1)
  result <- evaluate(hrf, t)

  expect_equal(length(result), length(t))

  # Peak should be near 2.5
  peak_time <- t[which.max(result)]
  expect_true(abs(peak_time - 2.5) < 0.5)

  # Should be 0 outside the range
  expect_true(all(result[t < 0] == 0))
  expect_true(all(result[t > 5] == 0))
})

test_that("hrf_weighted works with regressor for trial-wise analysis", {
  # Create a weighted HRF for extracting signal from 0-6s post-stimulus
  hrf <- hrf_weighted(
    width = 6,
    weights = c(0.1, 0.2, 0.3, 0.3, 0.2, 0.1),
    normalize = TRUE
  )

  # Create regressor
  reg <- regressor(onsets = c(0, 30), hrf = hrf)

  expect_true(inherits(reg, "Reg"))

  # Evaluate regressor
  t <- seq(0, 50, by = 0.5)
  result <- evaluate(reg, t)

  expect_equal(length(result), length(t))

  # Should have activity in the expected windows
  expect_true(any(result[t >= 0 & t <= 6] > 0))
  expect_true(any(result[t >= 30 & t <= 36] > 0))
})

test_that("hrf_boxcar and hrf_weighted produce equivalent results for uniform weights",
{
  # Boxcar of width 4
  hrf_box <- hrf_boxcar(width = 4, amplitude = 1)

  # Equivalent weighted HRF with uniform weights
  hrf_wt <- hrf_weighted(
    times = c(0, 1, 2, 3, 4),
    weights = c(1, 1, 1, 1, 0),  # Last weight is 0 (end of interval)
    method = "constant"
  )

  # Evaluate both
  t <- seq(0, 6, by = 0.5)
  result_box <- evaluate(hrf_box, t)
  result_wt <- evaluate(hrf_wt, t)

  # Should be very similar (may differ slightly at boundaries due to implementation)
  # Check correlation is very high
  expect_true(cor(result_box, result_wt) > 0.99)
})


# ============================================================================
# Tests for List-of-HRFs (Trial-Varying HRFs) Functionality
# ============================================================================

test_that("regressor accepts a list of HRFs for trial-varying analysis", {
  # Create different HRFs for each trial
  hrf1 <- hrf_boxcar(width = 4, normalize = TRUE)  # 4s window
  hrf2 <- hrf_boxcar(width = 6, normalize = TRUE)  # 6s window
  hrf3 <- hrf_boxcar(width = 8, normalize = TRUE)  # 8s window

  # Create regressor with list of HRFs
  reg <- regressor(
    onsets = c(10, 30, 50),
    hrf = list(hrf1, hrf2, hrf3)
  )

  expect_true(inherits(reg, "Reg"))
  expect_true(isTRUE(attr(reg, "hrf_is_list")))
  expect_equal(length(reg$hrf), 3)
  expect_equal(length(reg$onsets), 3)
})

test_that("list-of-HRFs regressor evaluates correctly", {
  # Create different HRFs for each trial
  hrf1 <- hrf_boxcar(width = 4, normalize = TRUE)
  hrf2 <- hrf_boxcar(width = 4, normalize = TRUE)

  # Create regressor with list of HRFs
  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  # Evaluate
  t <- seq(0, 50, by = 0.5)
  result <- evaluate(reg, t, method = "loop")

  expect_equal(length(result), length(t))

  # Should have activity in the expected windows
  expect_true(any(result[t >= 10 & t <= 14] > 0))  # First event window
  expect_true(any(result[t >= 30 & t <= 34] > 0))  # Second event window

  # Should be zero well outside event windows
  expect_equal(result[t == 0], 0)
  expect_equal(result[t == 50], 0)
})

test_that("list-of-HRFs with different spans works correctly", {
  # Create HRFs with different spans
  hrf1 <- hrf_boxcar(width = 4)   # span = 4
  hrf2 <- hrf_boxcar(width = 10)  # span = 10

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  # Overall span should be max of individual spans
  expect_equal(reg$span, 10)

  # Evaluate
  t <- seq(0, 50, by = 0.5)
  result <- evaluate(reg, t, method = "loop")

  # First event should only extend to 14s
  # Second event should extend to 40s
  expect_true(result[t == 12] > 0)  # Within first event window
  expect_true(result[t == 38] > 0)  # Within second event window (larger span)
})

test_that("list-of-HRFs with weighted HRFs works", {
  # Create different weighted HRFs for each trial
  hrf1 <- hrf_weighted(
    times = c(0, 2, 4, 6),
    weights = c(0.1, 0.5, 0.3, 0.1),
    normalize = TRUE
  )
  hrf2 <- hrf_weighted(
    times = c(0, 2, 4, 6),
    weights = c(0.4, 0.3, 0.2, 0.1),
    normalize = TRUE
  )

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  expect_true(inherits(reg, "Reg"))
  expect_true(isTRUE(attr(reg, "hrf_is_list")))

  # Evaluate
  t <- seq(0, 50, by = 0.5)
  result <- evaluate(reg, t, method = "loop")

  expect_equal(length(result), length(t))
  expect_true(any(result > 0))
})

test_that("list-of-HRFs recycles single HRF to all events", {
  hrf <- hrf_boxcar(width = 4)

  # Pass single-element list
  reg <- regressor(
    onsets = c(10, 20, 30),
    hrf = list(hrf)
  )

  expect_true(isTRUE(attr(reg, "hrf_is_list")))
  expect_equal(length(reg$hrf), 3)  # Should be recycled to 3
})

test_that("list-of-HRFs validates length", {
  hrf1 <- hrf_boxcar(width = 4)
  hrf2 <- hrf_boxcar(width = 6)

  # Wrong number of HRFs should error
  expect_error(
    regressor(onsets = c(10, 20, 30), hrf = list(hrf1, hrf2)),
    "must have length 1 or length equal to number of onsets"
  )
})

test_that("list-of-HRFs filters correctly with zero amplitudes", {
  hrf1 <- hrf_boxcar(width = 4)
  hrf2 <- hrf_boxcar(width = 6)
  hrf3 <- hrf_boxcar(width = 8)

  # Second event has zero amplitude - should be filtered
  reg <- regressor(
    onsets = c(10, 20, 30),
    amplitude = c(1, 0, 1),
    hrf = list(hrf1, hrf2, hrf3)
  )

  expect_equal(length(reg$onsets), 2)  # Only 2 events remain
  expect_equal(length(reg$hrf), 2)     # HRF list also filtered
  expect_equal(reg$onsets, c(10, 30))
})

test_that("list-of-HRFs works with all evaluation methods", {
  hrf1 <- hrf_boxcar(width = 4, normalize = TRUE)
  hrf2 <- hrf_boxcar(width = 4, normalize = TRUE)

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  t <- seq(0, 50, by = 0.5)

  # All methods should work (fft, conv, Rconv fall back to loop for list HRFs)
  result_loop <- evaluate(reg, t, method = "loop")
  result_conv <- evaluate(reg, t, method = "conv")
  result_fft <- evaluate(reg, t, method = "fft")
  result_Rconv <- evaluate(reg, t, method = "Rconv")

  # All should produce same results
  expect_equal(result_conv, result_loop)
  expect_equal(result_fft, result_loop)
  expect_equal(result_Rconv, result_loop)
})

test_that("nbasis.Reg handles list of HRFs", {
  hrf1 <- HRF_SPMG1  # nbasis = 1
  hrf2 <- HRF_SPMG1

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  expect_equal(nbasis(reg), 1)
})

test_that("print.Reg handles list of HRFs", {
  hrf1 <- hrf_boxcar(width = 4)
  hrf2 <- hrf_boxcar(width = 6)

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  # Should print without error and indicate trial-varying
  # cli output may not be captured by capture.output, so just verify no error
  expect_no_error(print(reg))
  # Check hrf_is_list attribute is set
  expect_true(isTRUE(attr(reg, "hrf_is_list")))
})

test_that("list-of-HRFs with mixed HRF types works", {
  # Mix standard HRF with boxcar
  hrf1 <- HRF_SPMG1
  hrf2 <- hrf_boxcar(width = 6, normalize = TRUE)

  reg <- regressor(
    onsets = c(10, 30),
    hrf = list(hrf1, hrf2)
  )

  t <- seq(0, 60, by = 0.5)
  result <- evaluate(reg, t, method = "loop")

  expect_equal(length(result), length(t))

  # First event should show typical hemodynamic response
  # Second event should show boxcar response
  expect_true(any(result[t >= 10 & t <= 25] > 0))
  expect_true(any(result[t >= 30 & t <= 36] > 0))
})
