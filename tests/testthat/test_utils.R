context("utility functions")

test_that("recycle_or_error recycles or errors correctly", {
  expect_equal(recycle_or_error(5, 3, "val"), rep(5,3))
  expect_equal(recycle_or_error(c(1,2,3), 3, "val"), c(1,2,3))
  expect_error(recycle_or_error(c(1,2), 3, "val"), "length 1 or 3")
})
