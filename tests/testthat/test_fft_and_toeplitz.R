context("FFT evaluation and HRF Toeplitz")

library(testthat)

# Helper HRF for Toeplitz and regressor tests
box_hrf_fun <- function(t) ifelse(t >= 0 & t <= 1, 1, 0)
BOX_HRF <- as_hrf(box_hrf_fun, name = "box", span = 1)


test_that("evaluate.Reg fft method matches conv", {
  reg <- regressor(onsets = c(0, 2), hrf = BOX_HRF, span = 1)
  grid <- seq(0, 4, by = 0.5)
  res_fft <- evaluate(reg, grid, method = "fft", precision = 0.1)
  res_conv <- evaluate(reg, grid, method = "conv", precision = 0.1)
  expect_equal(res_fft, res_conv, tolerance = 1e-6)
})

manual_toeplitz <- function(col, row) {
  nr <- length(col)
  nc <- length(row)
  out <- matrix(0, nr, nc)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (i >= j) {
        out[i, j] <- col[i - j + 1]
      } else {
        out[i, j] <- row[j - i + 1]
      }
    }
  }
  out
}


test_that("hrf_toeplitz constructs correct Toeplitz matrix", {
  time <- 0:2
  len <- 5
  H_dense <- hrf_toeplitz(BOX_HRF, time, len, sparse = FALSE)
  hreg <- BOX_HRF(time)
  col <- c(hreg, rep(0, len - length(hreg)))
  row <- c(hreg[1], rep(0, len - 1))
  expected <- manual_toeplitz(col, row)
  expect_true(inherits(H_dense, "Matrix"))
  expect_equal(as.matrix(H_dense), expected)

  H_sparse <- hrf_toeplitz(BOX_HRF, time, len, sparse = TRUE)
  expect_s4_class(H_sparse, "Matrix")
  expect_equal(as.matrix(H_sparse), expected)
})

