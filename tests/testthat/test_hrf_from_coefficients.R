context("hrf_from_coefficients")

test_that("hrf_from_coefficients combines basis correctly", {
  t <- seq(0, 20, by = 0.5)
  weights <- c(0.5, 2)
  h_combined <- hrf_from_coefficients(HRF_SPMG2, weights, name = "combined")
  expect_s3_class(h_combined, "HRF")
  expect_equal(nbasis(h_combined), 1L)
  expect_equal(attr(h_combined, "name"), "combined")

  expected <- as.numeric(HRF_SPMG2(t) %*% weights)
  result <- h_combined(t)
  expect_equal(result, expected)
})
