box_hrf_fun <- function(t) ifelse(t >= 0 & t <= 1, 1, 0)
BOX_HRF <- as_hrf(box_hrf_fun, name = "box", span = 1)

# Basic construction and evaluation with sampling frame

test_that("regressor_design produces expected design", {
  ons <- c(0, 1, 0, 2)
  fac <- factor(c("a", "b", "a", "b"))
  blk <- c(1L, 1L, 2L, 2L)
  sf  <- sampling_frame(blocklens = c(3, 3), TR = 1)
  dmat <- regressor_design(ons, fac, blk, sf, hrf = BOX_HRF, precision = 1)
  expect_true(is.matrix(dmat))
  expect_equal(nrow(dmat), length(samples(sf)))
  expect_equal(ncol(dmat), nbasis(BOX_HRF) * length(levels(fac)))
})

# Sparse output

test_that("regressor_design returns sparse matrix", {
  ons <- c(0, 2)
  fac <- factor(c("a", "b"))
  blk <- c(1L, 2L)
  sf  <- sampling_frame(blocklens = c(3, 3), TR = 1)
  dmat <- regressor_design(ons, fac, blk, sf, hrf = BOX_HRF,
                           precision = 1, sparse = TRUE)
  expect_true(inherits(dmat, "dgCMatrix"))
})

