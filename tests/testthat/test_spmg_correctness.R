test_that("canonical SPMG agrees with independent gamma densities", {
  t <- seq(0, 32, by = 0.001)
  expected <- dgamma(t, 6) - dgamma(t, 16) / 6
  expect_equal(hrf_spmg1(t), expected, tolerance = 1e-13)
  expect_equal(as.numeric(HRF_SPMG1(t)), expected, tolerance = 1e-13)
  expect_equal(-min(expected) / max(expected), 0.08891059, tolerance = 1e-7)
})

test_that("SPMG3 dispersion holds response mean and mass fixed", {
  t <- c(-1, 0, 0.1, 1, 3, 5, 7, 10, 16, 24, 32)
  base <- dgamma(t, 6) - dgamma(t, 16) / 6
  wider <- dgamma(t, 6 / 1.01, scale = 1.01) - dgamma(t, 16) / 6
  expected <- (base - wider) / 0.01
  actual <- HRF_SPMG3(t)[, 3]
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_gt(max(abs(actual - fmrihrf:::hrf_spmg1_second_deriv(t))), 0.02)
  expect_equal(integrate(function(x) HRF_SPMG3(x)[, 3], 0, 100)$value,
               0, tolerance = 1e-8)
  x <- seq(0.5, 24, by = 0.5)
  eps <- 1e-4
  numeric <- (HRF_SPMG3(x + eps) - HRF_SPMG3(x - eps)) / (2 * eps)
  expect_equal(deriv(HRF_SPMG3, x), numeric, tolerance = 1e-7)
})

test_that("SPMG fixes preserve fixed and legacy normalization policies", {
  h <- HRF_SPMG3
  t <- c(1, 3, 5, 9, 16, 24)
  raw <- h(t)
  expect_identical(normalize_hrf(h, "none"), h)
  ref <- seq(0, attr(h, "span"), length.out = round(attr(h, "span") * 50) + 1)
  rv <- h(ref)
  spm_ref <- seq(0, 32, length.out = 1600)
  factors <- list(spm = sum(h(spm_ref)[, 1]),
                  unit_peak = max(abs(rv[, 1])),
                  unit_integral = sum(diff(ref) * (head(rv[, 1], -1) + tail(rv[, 1], -1)) / 2),
                  unit_peak_per_basis = apply(abs(rv), 2, max))
  for (mode in names(factors)) {
    normalized <- normalize_hrf(h, mode)
    f <- factors[[mode]]
    expected <- if (length(f) == 1) raw / f else sweep(raw, 2, f, "/")
    expect_equal(normalized(t), expected, tolerance = 1e-12)
    expect_equal(normalized(t)[3, ], as.numeric(normalized(t[3])), tolerance = 1e-12)
    expect_equal(getHRF("spmg3", hrf_norm = mode)(t), expected, tolerance = 1e-12)
  }
  expect_equal(getHRF("spmg3", normalize = TRUE)(t), normalise_hrf(h)(t), tolerance = 1e-12)
})

test_that("loop block support agrees with an independent finite-kernel integral", {
  h <- as_hrf(function(t) as.numeric(t >= 0), name = "flat_tail", span = 2)
  reg <- regressor(0, hrf = h, duration = 2)
  grid <- c(0, 1, 2, 3, 4, 5)
  # Convolution of two unit-height boxes of width 2 is triangular.
  expected <- c(0, 1, 2, 1, 0, 0)
  for (method in c("loop", "conv")) {
    actual <- as.numeric(evaluate(reg, grid, precision = 0.001, method = method))
    expect_lt(max(abs(actual - expected)), 0.002)
  }
})
