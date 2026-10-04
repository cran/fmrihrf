# Regression tests for GitHub issue #44: fixed-scale HRF normalization.

.trapz44 <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

test_that("SPM normalization matches the Nilearn reference-grid convention", {
  h <- normalize_hrf(HRF_SPMG1, "spm")
  grid <- seq(0, 32, length.out = 1600L)

  expect_equal(sum(h(grid)), 1, tolerance = 1e-12)
  expect_equal(dim(h(6)), NULL)
})

test_that("fixed normalization is independent of the evaluation grid", {
  h <- normalize_hrf(HRF_SPMG1, "spm")
  factor1 <- HRF_SPMG1(6) / h(6)
  factor2 <- HRF_SPMG1(c(6, 8, 10))[1] / h(c(6, 8, 10))[1]

  expect_equal(factor1, factor2, tolerance = 1e-12)
})

test_that("unit peak and integral modes use the canonical basis scalar", {
  grid <- seq(0, attr(HRF_SPMG3, "span"), by = 0.02)
  raw <- HRF_SPMG3(grid)

  peak <- normalize_hrf(HRF_SPMG3, "unit_peak")(grid)
  integral <- normalize_hrf(HRF_SPMG3, "unit_integral")(grid)

  expect_equal(max(abs(peak[, 1L])), 1, tolerance = 1e-6)
  expect_equal(.trapz44(grid, integral[, 1L]), 1, tolerance = 1e-6)

  peak_factor <- raw[, 1L] / peak[, 1L]
  finite <- is.finite(peak_factor)
  expect_equal(peak[, 2L], raw[, 2L] / peak_factor[finite][1L],
               tolerance = 1e-12)
})

test_that("per-basis mode and normalise_hrf retain legacy semantics", {
  grid <- seq(0, attr(HRF_SPMG2, "span"), by = 0.02)
  explicit <- normalize_hrf(HRF_SPMG2, "unit_peak_per_basis")
  legacy <- normalise_hrf(HRF_SPMG2)

  expect_equal(legacy(grid), explicit(grid), tolerance = 0)
  expect_equal(apply(abs(explicit(grid)), 2L, max), c(1, 1),
               tolerance = 1e-6)
  expect_true(grepl("_norm", attr(legacy, "name"), fixed = TRUE))
})

test_that("construction entry points expose hrf_norm without changing defaults", {
  grid <- seq(0, 24, by = 0.1)

  expect_identical(gen_hrf(HRF_SPMG1)(grid), HRF_SPMG1(grid))
  expect_equal(sum(gen_hrf(HRF_SPMG1, hrf_norm = "spm")(
    seq(0, 32, length.out = 1600L)
  )), 1, tolerance = 1e-12)
  expect_equal(max(abs(getHRF("spmg1", hrf_norm = "unit_peak")(grid))),
               1, tolerance = 1e-3)

  expect_error(
    gen_hrf(HRF_SPMG1, normalize = TRUE, hrf_norm = "spm"),
    "either"
  )
  expect_error(
    getHRF("spmg1", normalize = TRUE, hrf_norm = "spm"),
    "either"
  )
})

test_that("normalization validates modes and unusable factors", {
  zero <- as_hrf(function(t) numeric(length(t)), name = "zero")

  expect_error(normalize_hrf(HRF_SPMG1), "mode")
  expect_error(normalize_hrf(HRF_SPMG1, "bogus"), "arg")
  expect_error(normalize_hrf(zero, "spm"), "not usable")
  expect_identical(normalize_hrf(HRF_SPMG1, "none"), HRF_SPMG1)
})
