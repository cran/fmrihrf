context("acquisition_onsets")

test_that("acquisition_onsets returns same values as samples(global=TRUE)", {
  # Single block
  sf1 <- sampling_frame(blocklens = 10, TR = 2)
  expect_equal(acquisition_onsets(sf1), samples(sf1, global = TRUE))
  
  # Multiple blocks
  sf2 <- sampling_frame(blocklens = c(10, 15, 20), TR = 2)
  expect_equal(acquisition_onsets(sf2), samples(sf2, global = TRUE))
  
  # Variable TR
  sf3 <- sampling_frame(blocklens = c(10, 10), TR = c(2, 1.5))
  expect_equal(acquisition_onsets(sf3), samples(sf3, global = TRUE))
})

test_that("acquisition_onsets handles default start_time correctly", {
  # Default start_time = TR/2
  sf <- sampling_frame(blocklens = 5, TR = 2)
  onsets <- acquisition_onsets(sf)
  
  # First onset should be at TR/2 = 1
  expect_equal(onsets[1], 1)
  
  # Spacing should be TR = 2
  expect_equal(diff(onsets), rep(2, 4))
})

test_that("acquisition_onsets handles custom start_time correctly", {
  # start_time = 0
  sf1 <- sampling_frame(blocklens = 5, TR = 2, start_time = 0)
  onsets1 <- acquisition_onsets(sf1)
  expect_equal(onsets1[1], 0)
  expect_equal(onsets1, c(0, 2, 4, 6, 8))
  
  # start_time = 3
  sf2 <- sampling_frame(blocklens = 5, TR = 2, start_time = 3)
  onsets2 <- acquisition_onsets(sf2)
  expect_equal(onsets2[1], 3)
  expect_equal(onsets2, c(3, 5, 7, 9, 11))
  
  # Variable start_time per block
  sf3 <- sampling_frame(blocklens = c(3, 3), TR = 2, start_time = c(0, 5))
  onsets3 <- acquisition_onsets(sf3)
  expect_equal(onsets3[1:3], c(0, 2, 4))      # Block 1
  expect_equal(onsets3[4:6], c(11, 13, 15))   # Block 2 starts at 6 + 5
})

test_that("acquisition_onsets handles variable TR correctly", {
  # Different TR per block
  sf <- sampling_frame(blocklens = c(5, 5), TR = c(2, 3))
  onsets <- acquisition_onsets(sf)
  
  # Block 1: TR=2, start_time=1
  expect_equal(onsets[1:5], c(1, 3, 5, 7, 9))
  
  # Block 2: TR=3, start_time=1.5, starts after block 1 (10 seconds)
  expect_equal(onsets[6:10], c(11.5, 14.5, 17.5, 20.5, 23.5))
})

test_that("acquisition_onsets works with multi-block designs", {
  # Three blocks
  sf <- sampling_frame(blocklens = c(100, 120, 80), TR = 2)
  onsets <- acquisition_onsets(sf)
  
  expect_equal(length(onsets), 300)  # Total scans
  
  # Check transitions between blocks
  # Block 1 ends at: 1 + 99*2 = 199
  expect_equal(onsets[100], 199)
  # Block 2 starts at: 200 + 1 = 201
  expect_equal(onsets[101], 201)
  
  # Block 2 ends at: 201 + 119*2 = 439
  expect_equal(onsets[220], 439)
  # Block 3 starts at: 440 + 1 = 441
  expect_equal(onsets[221], 441)
})

test_that("acquisition_onsets edge cases", {
  # Single scan
  sf1 <- sampling_frame(blocklens = 1, TR = 2)
  expect_equal(acquisition_onsets(sf1), 1)
  
  # Very short TR
  sf2 <- sampling_frame(blocklens = 5, TR = 0.5)
  onsets2 <- acquisition_onsets(sf2)
  expect_equal(onsets2, c(0.25, 0.75, 1.25, 1.75, 2.25))
  
  # Large number of blocks
  sf3 <- sampling_frame(blocklens = rep(10, 10), TR = 1)
  onsets3 <- acquisition_onsets(sf3)
  expect_equal(length(onsets3), 100)
})

test_that("acquisition_onsets matches expected timing for standard fMRI", {
  # Standard fMRI parameters
  sf <- sampling_frame(blocklens = c(150, 150), TR = 2, start_time = 0)
  onsets <- acquisition_onsets(sf)
  
  # Verify timing
  expect_equal(onsets[1], 0)      # First scan at t=0
  expect_equal(onsets[150], 298)  # Last scan of block 1
  expect_equal(onsets[151], 300)  # First scan of block 2
  expect_equal(onsets[300], 598)  # Last scan
  
  # Total duration
  expect_equal(max(onsets), (300-1) * 2)
})