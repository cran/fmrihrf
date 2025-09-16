context("penalty_matrix")

library(testthat)

test_that("bspline penalty is roughness-based", {
  R <- penalty_matrix(HRF_BSPLINE)
  expect_false(identical(R, diag(nbasis(HRF_BSPLINE))))
  expect_equal(nrow(R), nbasis(HRF_BSPLINE))
})

test_that("SPMG3 penalty shrinks derivatives", {
  R <- penalty_matrix(HRF_SPMG3, shrink_deriv = 4)
  expect_equal(R[1,1], 0)
  expect_equal(R[2,2], 4)
  expect_equal(R[3,3], 4)
})

test_that("fourier penalty increases with frequency", {
  fhrf <- hrf_fourier_generator(nbasis = 4, span = 24)
  R <- penalty_matrix(fhrf)
  expect_equal(diag(R), c(1,1,4,4))
})

test_that("default HRF penalty is identity", {
  R <- penalty_matrix(HRF_GAUSSIAN)
  expect_equal(R, diag(nbasis(HRF_GAUSSIAN)))
})
