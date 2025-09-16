context("Reg constructor and single_trial_regressor")

library(testthat)

# Verify Reg error handling

test_that("Reg validates inputs", {
  expect_error(Reg(onsets = -1), "onsets")
  expect_error(Reg(onsets = c(0, NA)), "onsets")
  expect_error(Reg(onsets = 1, span = 0), "span")
  expect_error(Reg(onsets = c(0,1), duration = c(1,2,3)), "duration")
  expect_error(Reg(onsets = c(0,1), amplitude = c(1,2,3)), "amplitude")
})

# Verify single_trial_regressor construction and errors

test_that("single_trial_regressor constructs single event", {
  st <- single_trial_regressor(onsets = 5, duration = 2, amplitude = 3)
  expect_s3_class(st, "Reg")
  expect_equal(length(st$onsets), 1)
  expect_equal(st$onsets, 5)
  expect_equal(st$duration, 2)
  expect_equal(st$amplitude, 3)
})

test_that("single_trial_regressor validates lengths", {
  expect_error(single_trial_regressor(onsets = c(1,2)), "onsets")
  expect_error(single_trial_regressor(onsets = 1, duration = c(1,2)), "duration")
  expect_error(single_trial_regressor(onsets = 1, amplitude = c(1,2)), "amplitude")
})

