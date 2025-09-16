context("hrf_lwu and related functions")

test_that("hrf_lwu computes response and normalisation", {
  t <- seq(0, 20, by = 0.5)
  raw <- hrf_lwu(t, tau = 6, sigma = 2, rho = 0.4)
  expect_equal(length(raw), length(t))
  expect_true(is.numeric(raw))

  # height normalisation
  norm <- hrf_lwu(t, tau = 6, sigma = 2, rho = 0.4, normalize = "height")
  expect_equal(max(abs(norm)), 1)

  # area normalisation currently behaves like none and warns
  expect_warning(area <- hrf_lwu(t, tau = 6, sigma = 2, rho = 0.4, normalize = "area"))
  expect_equal(area, raw)
})


test_that("hrf_basis_lwu returns derivatives", {
  theta <- c(tau = 6, sigma = 2, rho = 0.4)
  t <- seq(0, 20, by = 1)
  basis <- hrf_basis_lwu(theta, t)
  expect_equal(dim(basis), c(length(t), 4))
  expect_equal(basis[, "h0"], hrf_lwu(t, tau = theta["tau"], sigma = theta["sigma"], rho = theta["rho"]))

  # finite-difference check of derivative w.r.t tau at one point
  delta <- 1e-4
  t0 <- 4
  fd_tau <- (hrf_lwu(t0, tau = theta["tau"] + delta, sigma = theta["sigma"], rho = theta["rho"]) -
              hrf_lwu(t0, tau = theta["tau"] - delta, sigma = theta["sigma"], rho = theta["rho"])) / (2 * delta)
  idx <- which(t == t0)
  expect_equal(as.numeric(basis[idx, "d_tau"]), as.numeric(fd_tau), tolerance = 1e-3)
})
