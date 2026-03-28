context("Reg and regressor")

library(testthat)

# Helper HRF that returns 1 for t in [0,1] and 0 otherwise
box_hrf_fun <- function(t) {
  ifelse(t >= 0 & t <= 1, 1, 0)
}
BOX_HRF <- as_hrf(box_hrf_fun, name = "box", span = 1)


test_that("regressor constructs valid Reg objects", {
  reg <- regressor(onsets = c(0, 10, 20),
                   hrf = HRF_GAMMA,
                   duration = 2,
                   amplitude = c(1, 2, 3),
                   span = 30)

  expect_s3_class(reg, "Reg")
  expect_true(inherits(reg, "list"))
  expect_equal(reg$onsets, c(0, 10, 20))
  expect_equal(reg$duration, rep(2, 3))
  expect_equal(reg$amplitude, c(1, 2, 3))
  expect_identical(reg$hrf, HRF_GAMMA)
  # span should take from HRF when provided
  expect_equal(reg$span, attr(HRF_GAMMA, "span"))
  expect_false(attr(reg, "filtered_all"))
})


test_that("events with zero amplitude are filtered", {
  reg <- regressor(onsets = c(1, 2, 3), amplitude = c(1, 0, 2))
  expect_equal(reg$onsets, c(1, 3))
  expect_equal(reg$duration, c(0, 0))
  expect_equal(reg$amplitude, c(1, 2))
  expect_false(attr(reg, "filtered_all"))

  reg_empty <- regressor(onsets = c(1, 2), amplitude = c(0, 0))
  expect_length(reg_empty$onsets, 0)
  expect_true(attr(reg_empty, "filtered_all"))
})

test_that("NA inputs trigger errors", {
  # Single NA onset is treated as empty regressor (special case)
  expect_no_error(regressor(onsets = NA_real_))
  expect_error(regressor(onsets = c(1, NA)), "onsets")
  expect_error(regressor(onsets = 1, duration = NA_real_), "duration")
  expect_error(regressor(onsets = 1, amplitude = NA_real_), "amplitude")
  expect_error(regressor(onsets = 1, span = NA_real_), "span")
})


test_that("invalid inputs are rejected", {
  expect_error(regressor(onsets = c(-1, 1)), "onsets")
  expect_error(regressor(onsets = 1, duration = -2), "duration")
  expect_error(regressor(onsets = 1, span = 0), "span")
  expect_error(regressor(onsets = c(1, Inf)), "onsets")
  expect_error(regressor(onsets = 1, duration = Inf), "duration")
  expect_error(regressor(onsets = 1, amplitude = Inf), "amplitude")
  expect_error(regressor(onsets = 1, span = Inf), "span")
})


test_that("single_trial_regressor returns a length-1 Reg", {
  st <- single_trial_regressor(onsets = 5, duration = 2, amplitude = 3)
  expect_s3_class(st, "Reg")
  expect_equal(length(st$onsets), 1)
  expect_equal(st$onsets, 5)
  expect_equal(st$duration, 2)
  expect_equal(st$amplitude, 3)
})


test_that("shift.Reg shifts onsets correctly", {
  reg <- regressor(c(0, 2), hrf = HRF_SPMG1)
  shifted <- shift(reg, 5)
  expect_equal(shifted$onsets, c(5, 7))
  expect_identical(shifted$duration, reg$duration)
  expect_identical(shifted$amplitude, reg$amplitude)
})

test_that("shift.Reg accepts offset argument", {
  reg <- regressor(c(1, 3), hrf = HRF_SPMG1)
  shifted <- shift(reg, offset = 2)
  expect_equal(shifted$onsets, c(3, 5))
})

test_that("shift.Reg errors without shift specification", {
  reg <- regressor(0, hrf = HRF_SPMG1)
  expect_error(shift(reg), "shift_amount")
})


test_that("evaluate.Reg computes convolution correctly", {
  reg <- regressor(onsets = c(0, 2), hrf = BOX_HRF, span = 1)
  grid <- 0:4
  result <- evaluate(reg, grid, method = "conv", precision = 1)
  expect_equal(result, c(1, 1, 1, 1, 0))
})


test_that("Rconv handles non-zero constant durations", {
  reg <- regressor(onsets = c(0, 2), duration = 2, amplitude = c(1, 1),
                   hrf = BOX_HRF, span = 1)
  grid <- 0:6
  res_conv <- evaluate(reg, grid, method = "conv", precision = 1)
  res_rconv <- evaluate(reg, grid, method = "Rconv", precision = 1)
  expect_equal(res_rconv, res_conv)
})


test_that("unsorted grid triggers warning and sorted output", {
  reg <- regressor(onsets = 0, hrf = BOX_HRF, span = 1)
  expect_warning(out <- evaluate(reg, c(3, 0, 1), method = "conv", precision = 1))
  expect_equal(out, evaluate(reg, sort(c(3, 0, 1)), method = "conv", precision = 1))
})


test_that("evaluate.Reg validates grid and precision", {
  reg <- regressor(onsets = 0, hrf = BOX_HRF, span = 1)
  expect_error(evaluate(reg, numeric(0), method = "conv"), "grid")
  expect_error(evaluate(reg, c(0, NA), method = "conv"), "grid")
  expect_error(evaluate(reg, 0:1, precision = 0, method = "conv"), "precision")
  expect_error(evaluate(reg, 0:1, precision = -1, method = "conv"), "precision")
})

test_that("single-trial regressors with different durations normalize to peak 1", {
  grid <- seq(0, 40, by = 0.1)
  st_short <- single_trial_regressor(onsets = 10, duration = 1, hrf = HRF_SPMG1)
  st_long <- single_trial_regressor(onsets = 10, duration = 6, hrf = HRF_SPMG1)

  y_short <- evaluate(st_short, grid, method = "conv", precision = 0.1, normalize = TRUE)
  y_long <- evaluate(st_long, grid, method = "conv", precision = 0.1, normalize = TRUE)

  expect_equal(max(abs(y_short)), 1, tolerance = 1e-6)
  expect_equal(max(abs(y_long)), 1, tolerance = 1e-6)
})

test_that("evaluate.Reg normalize is consistent across methods", {
  reg <- regressor(onsets = c(2, 10), duration = c(1, 4), amplitude = c(1, 1), hrf = HRF_SPMG1)
  grid <- seq(0, 30, by = 0.25)

  y_conv <- evaluate(reg, grid, method = "conv", precision = 0.1, normalize = TRUE)
  y_loop <- evaluate(reg, grid, method = "loop", precision = 0.1, normalize = TRUE)

  expect_equal(max(abs(y_conv)), 1, tolerance = 1e-6)
  expect_equal(max(abs(y_loop)), 1, tolerance = 1e-6)
  expect_equal(y_conv, y_loop, tolerance = 0.05)
})

test_that("evaluate.Reg validates normalize argument", {
  reg <- regressor(onsets = 0, hrf = BOX_HRF, span = 1)
  expect_error(evaluate(reg, 0:1, normalize = NA), "normalize")
  expect_error(evaluate(reg, 0:1, normalize = 1), "normalize")
})
