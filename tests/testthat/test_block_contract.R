test_that("normalized blocked HRFs have query-independent scale", {
  h <- block_hrf(HRF_SPMG1, width = 4, normalize = TRUE)
  times <- c(1, 2.3, 7, 12, 25)
  whole <- h(times)
  # BLAS may use a different summation order for differently shaped queries.
  # Allow double-precision roundoff, not the old query-dependent rescaling.
  expect_equal(h(1), whole[1], tolerance = 64 * .Machine$double.eps)
  expect_equal(vapply(times, h, numeric(1)), whole, tolerance = 64 * .Machine$double.eps)
  expect_equal(c(h(times[1:2]), h(times[3:5])), whole, tolerance = 64 * .Machine$double.eps)
  expect_equal(h(rev(times)), rev(whole), tolerance = 64 * .Machine$double.eps)
  expect_equal(h(c(1, 1)), rep(whole[1], 2), tolerance = 64 * .Machine$double.eps)
  expect_lt(h(1), 0.01)
  expect_equal(attr(h, "span"), 28)
  expect_equal(attr(h, "params")$.normalize, TRUE)
  expect_equal(attr(h, "params")$A1, 1 / 120)

  # The compiled engine samples a full grid; the loop queries individual
  # event supports. Neither may change the normalization of the same object.
  reg <- regressor(0, h)
  for (method in c("conv", "loop")) {
    expect_equal(evaluate(reg, 1, precision = 0.01, method = method), h(1),
                 tolerance = 1e-12, info = method)
    expect_equal(evaluate(reg, c(1, 7), precision = 0.01, method = method),
                 h(c(1, 7)), tolerance = 1e-12, info = method)
  }
})

test_that("fixed block normalization preserves signs and scales bases independently", {
  triangle <- function(t) pmax(0, 1 - abs(t - 1))
  # Closed-form integral of the triangle, independent of package quadrature.
  primitive <- function(t) ifelse(t <= 0, 0,
    ifelse(t < 1, t^2 / 2, ifelse(t < 2, 1 - (2 - t)^2 / 2, 1)))
  basis <- as_hrf(function(t) cbind(
    signed = triangle(t) - 3 * triangle(t - 2),
    negative = -2 * triangle(t), zero = numeric(length(t))
  ), name = "signed triangles", nbasis = 3, span = 4)
  h <- block_hrf(basis, width = 4, precision = 0.02, normalize = TRUE)
  times <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8)
  signed_integral <- function(t) primitive(t) - 3 * primitive(t - 2)
  expected <- cbind(
    signed = (signed_integral(times) - signed_integral(times - 4)) / 3,
    negative = -(primitive(times) - primitive(times - 4)), zero = 0
  )
  # The signed column's absolute peak is at t = 6, beyond the original span
  # of 4. It is negative: normalization must neither truncate nor flip it.
  expect_equal(h(times), expected, tolerance = 1e-12)
  expect_equal(h(6)[1, ], c(signed = -1, negative = 0, zero = 0),
               tolerance = 1e-12)
  expect_equal(do.call(rbind, lapply(times, h)), h(times), tolerance = 64 * .Machine$double.eps)
  expect_equal(rbind(h(times[1:4]), h(times[5:10])), h(times), tolerance = 64 * .Machine$double.eps)
  expect_equal(dim(h(6)), c(1L, 3L))
  expect_equal(attr(h, "span"), 8)
  expect_true(all(is.finite(h(times))))
  expect_equal(h(times)[, "zero"], rep(0, length(times)))
})

test_that("every positive block width is integrated across the precision boundary", {
  flat <- as_hrf(function(t) as.numeric(t >= 0 & t <= 4), span = 4)
  linear <- as_hrf(function(t) ifelse(t >= 0 & t <= 4, t, 0), span = 4)
  widths <- c(1e-6, 0.05, 0.0999, 0.1, 0.1001)
  for (width in widths) {
    for (precision in c(0.1, 0.01)) {
      label <- paste("width", width, "precision", precision)
      expect_equal(block_hrf(flat, width, precision)(1), width,
                   tolerance = 1e-12, info = label)
      expect_equal(block_hrf(flat, width, precision, summate = FALSE)(1), 1,
                   tolerance = 1e-12, info = label)
      # Integral from 0 to width of (t-u) is t*width - width^2/2.
      times <- c(1, 2)
      expected <- times * width - width^2 / 2
      expect_equal(block_hrf(linear, width, precision)(times), expected,
                   tolerance = 1e-12, info = label)
      expect_equal(block_hrf(linear, width, precision, summate = FALSE)(times),
                   times - width / 2, tolerance = 1e-12, info = label)
    }
  }
  expect_equal(block_hrf(linear, 0, precision = 0.1)(c(1, 2)), c(1, 2))
  expect_equal(block_hrf(linear, 0, summate = FALSE)(c(1, 2)), c(1, 2))
  # At width = precision the model is continuous; only the quadrature changes.
  below <- block_hrf(flat, 0.1 - 1e-8, precision = 0.1)(1)
  above <- block_hrf(flat, 0.1 + 1e-8, precision = 0.1)(1)
  expect_equal(above - below, 2e-8, tolerance = 1e-12)
})

test_that("small blocks preserve multi-basis integration and decay", {
  h <- as_hrf(function(t) cbind(t, -2 * t, 0 * t), nbasis = 3, span = 4)
  width <- 0.05
  times <- c(1, 2)
  exact <- times * width - width^2 / 2
  expect_equal(unname(block_hrf(h, width, precision = 0.1)(times)),
               unname(cbind(exact, -2 * exact, 0 * exact)),
               tolerance = 1e-12)

  flat <- as_hrf(function(t) as.numeric(t >= 0 & t <= 4), span = 4)
  # Exponential decay has a closed-form integral, and both coarse and fine
  # trapezoid rules must approximate that integral rather than an impulse.
  half_life <- 0.7
  exact <- (1 - exp(-log(2) * width / half_life)) / (log(2) / half_life)
  coarse <- block_hrf(flat, width, precision = 0.1, half_life = half_life)(1)
  fine <- block_hrf(flat, width, precision = 0.001, half_life = half_life)(1)
  expect_lt(abs(coarse - exact), 1.1e-5)
  expect_lt(abs(fine - exact), 5e-9)
  # Preserve duration-averaging with attenuation; do not silently change to
  # renormalized exponential weights.
  expect_equal(block_hrf(flat, width, precision = 0.001,
                         half_life = half_life, summate = FALSE)(1),
               fine / width, tolerance = 1e-12)
})
