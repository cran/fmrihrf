test_that("deriv works for SPMG1 HRF", {
  t <- seq(0, 20, by = 0.5)
  
  # Test SPMG1 derivative
  deriv_vals <- deriv(HRF_SPMG1, t)
  
  # Should return a numeric vector
  expect_true(is.numeric(deriv_vals))
  expect_equal(length(deriv_vals), length(t))
  
  # Derivative should be 0 at t=0
  expect_equal(deriv_vals[1], 0)
  
  # Check that it matches the analytic derivative
  params <- attr(HRF_SPMG1, "params")
  if (is.null(params)) {
    params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  }
  expected <- hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  expect_equal(deriv_vals, expected)
})

test_that("deriv works for SPMG2 HRF", {
  t <- seq(0, 20, by = 0.5)
  
  # Test SPMG2 derivative
  deriv_vals <- deriv(HRF_SPMG2, t)
  
  # Should return a matrix with 2 columns
  expect_true(is.matrix(deriv_vals))
  expect_equal(ncol(deriv_vals), 2)
  expect_equal(nrow(deriv_vals), length(t))
  
  # First column should match SPMG1 derivative
  params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  expected_col1 <- hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  expect_equal(deriv_vals[, 1], expected_col1)
  
  # Second column should match second derivative
  expected_col2 <- hrf_spmg1_second_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  expect_equal(deriv_vals[, 2], expected_col2)
})

test_that("deriv works for SPMG3 HRF", {
  t <- seq(0, 20, by = 0.5)
  
  # Test SPMG3 derivative
  deriv_vals <- deriv(HRF_SPMG3, t)
  
  # Should return a matrix with 3 columns
  expect_true(is.matrix(deriv_vals))
  expect_equal(ncol(deriv_vals), 3)
  expect_equal(nrow(deriv_vals), length(t))
  
  # First two columns should match SPMG2 derivatives
  params <- list(P1 = 5, P2 = 15, A1 = 0.0833)
  expected_col1 <- hrf_spmg1_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  expect_equal(deriv_vals[, 1], expected_col1)
  
  expected_col2 <- hrf_spmg1_second_deriv(t, P1 = params$P1, P2 = params$P2, A1 = params$A1)
  expect_equal(deriv_vals[, 2], expected_col2)
  
  # Third column should be numeric (third derivative)
  expect_true(all(is.numeric(deriv_vals[, 3])))
})

test_that("deriv works with numeric differentiation for other HRF types", {
  skip_if_not_installed("numDeriv")
  
  t <- seq(0, 10, by = 0.5)
  
  # Test Gaussian HRF
  deriv_gauss <- deriv(HRF_GAUSSIAN, t)
  expect_true(is.numeric(deriv_gauss))
  expect_equal(length(deriv_gauss), length(t))
  
  # Test Gamma HRF
  deriv_gamma <- deriv(HRF_GAMMA, t)
  expect_true(is.numeric(deriv_gamma))
  expect_equal(length(deriv_gamma), length(t))
})

test_that("deriv works for multi-basis HRFs", {
  skip_if_not_installed("numDeriv")
  
  t <- seq(0, 20, by = 1)
  
  # Test B-spline HRF
  deriv_bspline <- deriv(HRF_BSPLINE, t)
  expect_true(is.matrix(deriv_bspline))
  expect_equal(ncol(deriv_bspline), nbasis(HRF_BSPLINE))
  expect_equal(nrow(deriv_bspline), length(t))
})

test_that("deriv handles edge cases", {
  # Test with single time point
  deriv_single <- deriv(HRF_SPMG1, 5)
  expect_true(is.numeric(deriv_single))
  expect_equal(length(deriv_single), 1)
  
  # Test with empty time vector
  deriv_empty <- deriv(HRF_SPMG1, numeric(0))
  expect_true(is.numeric(deriv_empty))
  expect_equal(length(deriv_empty), 0)
  
  # Test with negative times
  t_neg <- c(-5, -2, 0, 2, 5)
  deriv_neg <- deriv(HRF_SPMG1, t_neg)
  expect_true(is.numeric(deriv_neg))
  expect_equal(length(deriv_neg), length(t_neg))
  # Derivative should be 0 for negative times
  expect_equal(deriv_neg[1:2], c(0, 0))
})