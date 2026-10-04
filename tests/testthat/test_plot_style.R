test_that("hrf_palette returns the documented palettes", {
  expect_length(hrf_palette(), 6)
  expect_identical(hrf_palette(3), hrf_palette()[1:3])
  expect_length(hrf_palette(12, "ordered"), 12)
  expect_identical(hrf_palette(0), character(0))
  # More than six categorical colours fall back to the ordered ramp
  expect_identical(hrf_palette(8), hrf_palette(8, "ordered"))
  # Brick red is reserved for the categorical palette
  expect_false(hrf_palette(1) %in% hrf_palette(12, "ordered"))
  expect_error(hrf_palette(-1), "non-negative")
})

test_that("palette colours meet 3:1 contrast on light and dark surfaces", {
  lum <- function(h) {
    m <- grDevices::col2rgb(h) / 255
    m <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
    0.2126 * m[1, ] + 0.7152 * m[2, ] + 0.0722 * m[3, ]
  }
  ratio <- function(a, b) (pmax(lum(a), lum(b)) + 0.05) / (pmin(lum(a), lum(b)) + 0.05)
  for (cols in list(hrf_palette(), hrf_palette(12, "ordered"))) {
    expect_true(all(ratio(cols, "#ffffff") >= 3))
    expect_true(all(ratio(cols, "#14171c") >= 3))
  }
})

test_that("ggplot2 scales use hrf_palette", {
  skip_if_not_installed("ggplot2")
  sc <- scale_colour_hrf()
  expect_s3_class(sc, "ScaleDiscrete")
  expect_identical(sc$palette(3), hrf_palette(3))
  expect_identical(scale_fill_hrf("ordered")$palette(5), hrf_palette(5, "ordered"))
})

test_that("plot_hrfs expands basis sets and draws a reference curve", {
  skip_if_not_installed("ggplot2")
  time <- seq(0, 20, by = 0.5)
  d <- plot_hrfs(HRF_SPMG3, time = time, draw = FALSE)
  expect_equal(levels(d$HRF), c("B1", "B2", "B3"))
  expect_equal(d$response, as.vector(HRF_SPMG3(time)))
  expect_equal(nlevels(plot_hrfs(HRF_SPMG3, time = time, basis = "first", draw = FALSE)$HRF), 1)

  d <- plot_hrfs(HRF_GAMMA, time = time, reference = HRF_SPMG1, draw = FALSE)
  p <- attr(d, "plot")
  ref_layer <- Filter(function(l) identical(l$aes_params$linetype, "22"), p$layers)
  expect_length(ref_layer, 1)
  expect_equal(ref_layer[[1]]$data$response, HRF_SPMG1(time))
  expect_error(plot_hrfs(HRF_GAMMA, reference = 1), "HRF object")
})

test_that("plot_regressors handles RegSets, stacked onsets and precision", {
  skip_if_not_installed("ggplot2")
  set <- regressor_set(c(2, 10, 20, 30), factor(c("a", "b", "a", "b")))
  grid <- seq(0, 50, by = 0.1)
  d <- plot_regressors(set, grid = grid, layout = "stack", draw = FALSE)
  expect_equal(levels(d$Regressor), c("a", "b"))
  # In a stacked layout each panel is marked with its own events
  rects <- Filter(function(l) inherits(l$geom, "GeomRect"), attr(d, "plot")$layers)[[1]]$data
  expect_equal(sort(rects$onset[rects$series == "a"]), c(2, 20))
  expect_equal(sort(rects$onset[rects$series == "b"]), c(10, 30))

  # Onsets outside the grid are not drawn
  r <- regressor(c(5, 100), HRF_SPMG1)
  d <- plot_regressors(r, grid = seq(0, 40, by = 0.5), show_onsets = TRUE, draw = FALSE)
  rects <- Filter(function(l) inherits(l$geom, "GeomRect"), attr(d, "plot")$layers)[[1]]$data
  expect_equal(rects$onset, 5)

  # Default precision follows the grid spacing
  box <- regressor(10, hrf_boxcar(width = 4))
  fine <- seq(0, 20, by = 0.05)
  d <- plot_regressors(box, grid = fine, draw = FALSE)
  expect_equal(d$response, evaluate(box, fine, precision = 0.05))
})

test_that("results returned inside knitr subset to plain data frames", {
  skip_if_not_installed("ggplot2")
  old <- options(knitr.in.progress = TRUE)
  on.exit(options(old), add = TRUE)
  d <- plot_hrfs(HRF_SPMG1, HRF_GAMMA, time = seq(0, 10, by = 1))
  expect_s3_class(d, "fmrihrf_plot")
  h <- head(d)
  expect_false(inherits(h, "fmrihrf_plot"))
  expect_null(attr(h, "plot"))
  expect_equal(nrow(h), 6)
})

test_that("base plot methods use the shared palette and stack basis regressors", {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- plot(regressor(c(5, 30), HRF_SPMG3), grid = seq(0, 50, by = 0.5))
  expect_equal(names(d), c("time", "basis_1", "basis_2", "basis_3"))
  d <- plot(HRF_SPMG2)
  expect_equal(ncol(d), 3)
})

test_that("show_onsets = 'first' marks only the first regressor, also when stacked", {
  skip_if_not_installed("ggplot2")
  r1 <- regressor(c(2, 20), HRF_SPMG1)
  r2 <- regressor(c(10, 30), HRF_SPMG1)
  d <- plot_regressors(r1, r2, labels = c("one", "two"), grid = seq(0, 50, by = 0.1),
                       layout = "stack", show_onsets = "first", draw = FALSE)
  rects <- Filter(function(l) inherits(l$geom, "GeomRect"), attr(d, "plot")$layers)[[1]]$data
  expect_equal(unique(as.character(rects$series)), "one")
  expect_equal(sort(rects$onset), c(2, 20))
})

test_that("reference curves follow normalize", {
  skip_if_not_installed("ggplot2")
  time <- seq(0, 20, by = 0.1)
  d <- plot_hrfs(HRF_GAMMA, time = time, normalize = TRUE, reference = HRF_SPMG1,
                 draw = FALSE)
  ref <- Filter(function(l) identical(l$aes_params$linetype, "22"), attr(d, "plot")$layers)[[1]]$data
  expect_equal(max(ref$response), 1)
  expect_equal(max(d$response), 1)
})

test_that("modifying a knitr result returns a plain data frame", {
  skip_if_not_installed("ggplot2")
  old <- options(knitr.in.progress = TRUE)
  on.exit(options(old), add = TRUE)
  d <- plot_hrfs(HRF_SPMG1, time = seq(0, 10, by = 1))
  d2 <- d
  d2$extra <- 1
  expect_false(inherits(d2, "fmrihrf_plot"))
  expect_null(attr(d2, "plot"))
  d3 <- d
  d3[["response"]] <- 0
  expect_false(inherits(d3, "fmrihrf_plot"))
  expect_false(inherits(as.data.frame(d), "fmrihrf_plot"))
  expect_true(inherits(d, "fmrihrf_plot"))
})

test_that("y breaks are evenly spaced, within budget, and clear of event marks", {
  br <- fmrihrf:::.budget_breaks
  ranges <- list(c(-0.09, 0.16), c(-0.26, 0.6), c(-1, 1), c(-0.3, 0.65),
                 c(-2.4, 3.5), c(-0.42, 1.55), c(0, 1), c(-0.02, 0.175), c(0, 0.2))
  for (lim in ranges) {
    for (stack in c(FALSE, TRUE)) {
      b <- br(lim, stack = stack)
      expect_lte(length(b), if (stack) 3 else 5)
      if (!stack && length(b) == 5) expect_true(any(b < 0) && any(b > 0))
      if (length(b) > 2) expect_equal(diff(range(diff(b))), 0, tolerance = 1e-9)
      if (length(b) > 1) expect_gte(min(diff(b)), (if (stack) 0.22 else 0.12) * diff(lim) - 1e-9)
    }
  }
  # Both signs labelled where an even grid allows it
  for (lim in list(c(-1, 1), c(-0.26, 0.6), c(-2.4, 3.5), c(-0.09, 0.16))) {
    b <- br(lim)
    expect_true(any(b < 0) && any(b > 0) && any(b == 0))
  }
  expect_equal(br(c(0, 1)), c(0, 0.5, 1))
  # No labels inside the band reserved for event marks
  expect_true(all(br(c(-0.06, 0.175), floor = 0.2) >= 0))
})

test_that("stacked panels keep tick labels out of the event-mark band", {
  skip_if_not_installed("ggplot2")
  # Amplitudes -2.5 and -0.9 put a label inside the mark band with the
  # earlier 14% floor for stacked panels.
  r1 <- regressor(c(0, 12, 24), HRF_SPMG1, amplitude = c(-2.5, 1, 2))
  r2 <- regressor(c(0, 12, 24), HRF_SPMG1, amplitude = c(-0.9, 1, 2))
  d <- plot_regressors(r1, r2, grid = seq(0, 50, by = 0.1), layout = "stack",
                       draw = FALSE)
  b <- ggplot2::ggplot_build(attr(d, "plot"))
  rects <- b$data[[which(vapply(attr(d, "plot")$layers,
                                function(l) inherits(l$geom, "GeomRect"), logical(1)))]]
  for (panel in unique(rects$PANEL)) {
    top_of_marks <- max(rects$ymax[rects$PANEL == panel])
    breaks <- b$layout$panel_params[[as.integer(panel)]]$y$breaks
    breaks <- breaks[!is.na(breaks)]
    expect_true(all(breaks > top_of_marks))
  }
})


test_that("stacked mark band holds for adversarial and fixed-scale panels", {
  skip_if_not_installed("ggplot2")
  band_ok <- function(d) {
    p <- attr(d, "plot")
    b <- ggplot2::ggplot_build(p)
    i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1)))
    rc <- b$data[[i]]
    all(vapply(unique(rc$PANEL), function(panel) {
      top <- max(rc$ymax[rc$PANEL == panel])
      br <- b$layout$panel_params[[as.integer(panel)]]$y$breaks
      all(br[!is.na(br)] > top)
    }, logical(1)))
  }
  # Minimal case from review: -0.2 fell inside marks spanning [-0.238, -0.199]
  r <- regressor(c(1, 15, 29), HRF_SPMG1, amplitude = c(0.34, -1.03, -0.28))
  expect_true(band_ok(plot_regressors(r, grid = seq(0, 60, by = 0.1),
                                      layout = "stack", draw = FALSE)))
  set.seed(1)
  for (k in 1:25) {
    a <- round(stats::runif(3, -1.2, 1.2), 2)
    r1 <- regressor(c(1, 15, 29), HRF_SPMG1, amplitude = a)
    r2 <- regressor(c(5, 20), HRF_SPMG1, amplitude = rev(a)[1:2])
    for (sc in c("free_y", "fixed")) {
      d <- plot_regressors(r1, r2, grid = seq(0, 60, by = 0.1), layout = "stack",
                           scales = sc, draw = FALSE)
      expect_true(band_ok(d), info = paste(sc, paste(a, collapse = ",")))
    }
  }
})
