test_that("native method errors return safely through the Rcpp wrapper", {
  for (i in seq_len(3)) {
    expect_error(
      fmrihrf:::evaluate_regressor_cpp(
        grid = c(0, 1), onsets = 0, durations = 0, amplitudes = 1,
        hrf_matrix = matrix(c(0, 1, 0), ncol = 1), hrf_span = 2,
        precision = 1, method = "unsupported"
      ),
      "only 'conv' is supported", fixed = TRUE
    )
  }
  # An error must leave the R session and the native evaluator usable.
  h <- as_hrf(function(t) pmax(0, 1 - abs(t - 1)), span = 2)
  reg <- regressor(0, h)
  expect_equal(evaluate(reg, c(0, 1, 2), precision = 1), c(0, 1, 0))
})
