library(testthat)

# simple generator that records scale
make_gauss <- function(mean=0, sd=1, scale=1) {
  as_hrf(function(t) scale * dnorm(t, mean, sd),
         name = paste0("gauss_", mean),
         nbasis = 1L)
}


test_that("hrf_library forwards extra arguments", {
  pgrid <- data.frame(mean = c(0, 2), sd = c(1, 1))
  lib <- hrf_library(make_gauss, pgrid, scale = 2)
  expect_s3_class(lib, "HRF")
  expect_equal(nbasis(lib), 2)
  t <- c(0, 1)
  expected <- cbind(2 * dnorm(t, 0, 1), 2 * dnorm(t, 2, 1))
  expect_equal(lib(t), expected)
})


test_that("duplicate parameter names cause an error", {
  pgrid <- data.frame(mean = 0, sd = 1, scale = 1)
  expect_error(hrf_library(make_gauss, pgrid, scale = 2))
})
