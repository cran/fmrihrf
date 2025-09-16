context("Sampling Frame")

test_that("sampling_frame constructor works correctly", {
  # Basic construction
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  expect_s3_class(sframe, "sampling_frame")
  expect_equal(length(sframe$blocklens), 2)
  expect_equal(sframe$TR, c(2, 2))
  expect_equal(sframe$start_time, c(1, 1))
  
  # Test with different TRs per block
  sframe2 <- sampling_frame(blocklens = c(100, 200), TR = c(2, 1.5))
  expect_equal(sframe2$TR, c(2, 1.5))
  
  # Test input validation
  expect_error(sampling_frame(blocklens = c(-1, 100), TR = 2),
              "Block lengths must be positive")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = -1),
              "TR .* must be positive")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = 2, precision = 3),
              "Precision must be positive and less than")

  expect_error(sampling_frame(blocklens = c("a", 100), TR = 2),
               "numeric")
  expect_error(sampling_frame(blocklens = c(100, NA), TR = 2),
               "non-NA")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = c(2, NA)),
               "non-NA")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = "a"),
               "numeric")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = 2,
                              start_time = c(0, NA)),
               "non-NA")

  expect_error(sampling_frame(blocklens = c(100, 100), TR = c(2, 2, 2)),
              "TR must have length 1 or match the number of blocks")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = 2, start_time = c(0, 0, 0)),
              "start_time must have length 1 or match the number of blocks")

})

test_that("samples.sampling_frame works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Test relative timing
  rel_samples <- samples(sframe, global = FALSE)
  expect_equal(length(rel_samples), 200)
  expect_equal(rel_samples[1:5], c(1, 3, 5, 7, 9))
  
  # Test global timing
  glob_samples <- samples(sframe, global = TRUE)
  expect_equal(length(glob_samples), 200)
  expect_equal(glob_samples[101] - glob_samples[100], 2)  # Check TR spacing
  
  # Test block selection
  block1_samples <- samples(sframe, blockids = 1)
  expect_equal(length(block1_samples), 100)
  
  # Test memoization
  samples2 <- samples(sframe, global = FALSE)
  expect_identical(rel_samples, samples2)  # Should return cached result
})

test_that("global_onsets works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Test basic functionality
  onsets <- c(10, 20)
  blockids <- c(1, 2)
  global_times <- global_onsets(sframe, onsets, blockids)
  expect_equal(length(global_times), 2)
  expect_equal(global_times[1], 10)  # First block onset unchanged
  expect_equal(global_times[2], 220)  # Second block onset = 200 (block1 duration) + 20
  
  # Test error conditions for non-integer block ids
  expect_error(global_onsets(sframe, onsets, c(1.5, 2)),
               "blockids must be whole numbers")
  expect_error(global_onsets(sframe, onsets, c(1, NA)),
               "blockids must be whole numbers")
})

test_that("print.sampling_frame works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  expect_output(print(sframe), "Sampling Frame")
  expect_output(print(sframe), "Structure")
  expect_output(print(sframe), "Timing")
  expect_output(print(sframe), "Duration")
})

test_that("sampling_frame handles edge cases", {
  # Single block
  single_block <- sampling_frame(blocklens = 100, TR = 2)
  expect_equal(length(single_block$blocklens), 1)
  expect_equal(length(samples(single_block)), 100)
  
  # Very short block
  short_block <- sampling_frame(blocklens = c(1, 1), TR = 2)
  expect_equal(length(samples(short_block)), 2)
  
  # Different start times
  custom_starts <- sampling_frame(blocklens = c(100, 100), 
                                TR = 2, 
                                start_time = c(0, 5))
  expect_equal(custom_starts$start_time, c(0, 5))
  
  # High precision
  high_prec <- sampling_frame(blocklens = c(10, 10), 
                            TR = 2, 
                            precision = 0.01)
  expect_equal(high_prec$precision, 0.01)
})

test_that("sampling_frame maintains temporal consistency", {
  sframe <- sampling_frame(blocklens = c(100, 100, 100), TR = 2)
  glob_samples <- samples(sframe, global = TRUE)
  
  # Check uniform spacing within blocks
  for (block in 1:3) {
    block_idx <- which(blockids(sframe) == block)
    diffs <- diff(glob_samples[block_idx])
    expect_true(all(abs(diffs - 2) < 1e-10))
  }
  
  # Check block transitions
  block_ends <- cumsum(sframe$blocklens)
  for (i in 1:(length(block_ends)-1)) {
    time_diff <- glob_samples[block_ends[i] + 1] - glob_samples[block_ends[i]]
    expect_equal(time_diff, 2)
  }
})

test_that("samples with blockids parameter works correctly for local and global timing", {
  sframe <- sampling_frame(blocklens = c(10, 10, 10), TR = 2)
  
  # Test 1: Verify blockids filtering for local timing
  block1_local <- samples(sframe, blockids = 1, global = FALSE)
  block2_local <- samples(sframe, blockids = 2, global = FALSE)
  block3_local <- samples(sframe, blockids = 3, global = FALSE)
  
  expect_equal(length(block1_local), 10)
  expect_equal(length(block2_local), 10)
  expect_equal(length(block3_local), 10)
  
  # All blocks should have same local timing pattern
  expect_equal(block1_local, block2_local)
  expect_equal(block2_local, block3_local)
  
  # Test 2: Verify blockids filtering for global timing
  block1_global <- samples(sframe, blockids = 1, global = TRUE)
  block2_global <- samples(sframe, blockids = 2, global = TRUE)
  block3_global <- samples(sframe, blockids = 3, global = TRUE)
  
  expect_equal(length(block1_global), 10)
  expect_equal(length(block2_global), 10)
  expect_equal(length(block3_global), 10)
  
  # Verify correct global offsets
  expect_equal(block1_global[1], 1)      # First sample of block 1
  expect_equal(block2_global[1], 21)     # First sample of block 2 (after 10*2=20s)
  expect_equal(block3_global[1], 41)     # First sample of block 3 (after 20*2=40s)
  
  # Test 3: Multiple blocks selection
  blocks12_global <- samples(sframe, blockids = c(1, 2), global = TRUE)
  expect_equal(length(blocks12_global), 20)
  expect_equal(blocks12_global[1:10], block1_global)
  expect_equal(blocks12_global[11:20], block2_global)
  
  # Test 4: Non-sequential block selection
  blocks13_global <- samples(sframe, blockids = c(1, 3), global = TRUE)
  expect_equal(length(blocks13_global), 20)
  expect_equal(blocks13_global[1:10], block1_global)
  expect_equal(blocks13_global[11:20], block3_global)
  
  # Test 5: Ensure it does NOT return all samples when specific blocks requested
  all_samples <- samples(sframe, global = TRUE)
  expect_equal(length(all_samples), 30)
  expect_true(length(block1_global) < length(all_samples))
})

test_that("samples with varying TR and start_time per block work correctly", {
  # Different TR per block
  sframe_varTR <- sampling_frame(
    blocklens = c(10, 10, 10), 
    TR = c(2, 1.5, 3),
    start_time = c(1, 0.75, 1.5)
  )
  
  # Test local timing respects per-block TR
  block1_local <- samples(sframe_varTR, blockids = 1, global = FALSE)
  block2_local <- samples(sframe_varTR, blockids = 2, global = FALSE)
  block3_local <- samples(sframe_varTR, blockids = 3, global = FALSE)
  
  # Check spacing
  expect_equal(diff(block1_local)[1], 2)    # TR = 2
  expect_equal(diff(block2_local)[1], 1.5)  # TR = 1.5
  expect_equal(diff(block3_local)[1], 3)    # TR = 3
  
  # Check start times
  expect_equal(block1_local[1], 1)      # start_time = 1
  expect_equal(block2_local[1], 0.75)   # start_time = 0.75
  expect_equal(block3_local[1], 1.5)    # start_time = 1.5
  
  # Test global timing
  block1_global <- samples(sframe_varTR, blockids = 1, global = TRUE)
  block2_global <- samples(sframe_varTR, blockids = 2, global = TRUE)
  block3_global <- samples(sframe_varTR, blockids = 3, global = TRUE)
  
  # Block 1: starts at 0 + 1 = 1
  expect_equal(block1_global[1], 1)
  
  # Block 2: starts at 20 (10*2) + 0.75 = 20.75
  expect_equal(block2_global[1], 20.75)
  
  # Block 3: starts at 20 + 15 (10*1.5) + 1.5 = 36.5
  expect_equal(block3_global[1], 36.5)
})

test_that("edge cases for blockids parameter", {
  sframe <- sampling_frame(blocklens = c(5, 10, 15), TR = 1)
  
  # Empty blockids
  expect_equal(length(samples(sframe, blockids = integer(0))), 0)
  
  # Out of order blockids
  out_order <- samples(sframe, blockids = c(3, 1, 2), global = TRUE)
  expect_equal(length(out_order), 30)
  
  # Repeated blockids
  repeated <- samples(sframe, blockids = c(1, 1), global = TRUE)
  expect_equal(length(repeated), 10)  # Block 1 appears twice, so 5*2 = 10
  
  # Single block from middle
  middle_block <- samples(sframe, blockids = 2, global = TRUE)
  expect_equal(length(middle_block), 10)
  expect_equal(middle_block[1], 5.5)  # After block 1 (5 scans * 1s TR) + start_time
})