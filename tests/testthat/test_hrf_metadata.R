test_that("closed HRF constructors retain parameters without argument warnings", {
  expect_no_warning(box <- hrf_boxcar(4, amplitude = 2))
  expect_equal(attr(box, "params")$width, 4)
  expect_equal(box(c(-1, 0, 3, 4)), c(0, 2, 2, 0))

  expect_no_warning(weighted <- hrf_weighted(c(1, 2, 3), width = 2))
  expect_equal(attr(weighted, "params")$weights, c(1, 2, 3))

  expect_no_warning(lagged <- lag_hrf(HRF_SPMG1, 2))
  expect_equal(attr(lagged, "params")$A1, 1 / 120)
  expect_equal(attr(lagged, "params")$.lag, 2)
  expect_equal(lagged(c(2, 4, 8)), HRF_SPMG1(c(0, 2, 6)))

  expect_no_warning(blocked <- block_hrf(HRF_SPMG1, 2))
  expect_equal(attr(blocked, "params")$A1, 1 / 120)
  expect_equal(attr(blocked, "params")$.width, 2)

  expect_no_warning(combined <- hrf_from_coefficients(HRF_SPMG3, c(1, 2, 3)))
  expect_equal(attr(combined, "params")$coefficients, c(1, 2, 3))
  expect_equal(combined(c(1, 4, 8)), drop(HRF_SPMG3(c(1, 4, 8)) %*% c(1, 2, 3)))

  for (generator in list(hrf_tent_generator, hrf_fourier_generator,
                         hrf_daguerre_generator)) {
    expect_no_warning(basis <- generator(nbasis = 4))
    expect_true(length(attr(basis, "params")) > 0)
    expect_identical(attr(basis, "param_names"), names(attr(basis, "params")))
    expect_equal(dim(basis(c(1, 2, 3))), c(3, 4))
  }
})
