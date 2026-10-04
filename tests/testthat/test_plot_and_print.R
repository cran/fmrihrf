.with_pdf_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  value <- force(code)
  stopifnot(file.exists(path))
  value
}

test_that("plot.Reg returns the evaluated grid data", {
  reg <- regressor(onsets = c(2, 8), hrf = HRF_SPMG1)
  grid <- seq(0, 12, by = 0.5)

  plotted <- .with_pdf_device(
    plot(reg, grid = grid, show_onsets = FALSE)
  )

  expect_s3_class(plotted, "data.frame")
  expect_equal(names(plotted), c("time", "response"))
  expect_equal(plotted$time, grid)
  expect_equal(plotted$response, evaluate(reg, grid))
})

test_that("plot_regressors returns labelled long-form data", {
  reg1 <- regressor(onsets = c(2, 8), hrf = HRF_SPMG1)
  reg2 <- regressor(onsets = c(4, 10), hrf = HRF_GAMMA)
  grid <- seq(0, 12, by = 1)

  plotted <- .with_pdf_device(
    plot_regressors(reg1, reg2, grid = grid, labels = c("spm", "gamma"),
                    show_onsets = FALSE, use_ggplot = FALSE)
  )

  expect_s3_class(plotted, "data.frame")
  expect_equal(names(plotted), c("time", "Regressor", "response"))
  expect_equal(nrow(plotted), 2 * length(grid))
  expect_equal(levels(plotted$Regressor), c("spm", "gamma"))
})

test_that("plot_hrfs returns labelled long-form data", {
  time <- seq(0, 12, by = 1)

  plotted <- .with_pdf_device(
    plot_hrfs(HRF_SPMG1, HRF_GAMMA, time = time,
              labels = c("spm", "gamma"), use_ggplot = FALSE)
  )

  expect_s3_class(plotted, "data.frame")
  expect_equal(names(plotted), c("time", "HRF", "response"))
  expect_equal(nrow(plotted), 2 * length(time))
  expect_equal(levels(plotted$HRF), c("spm", "gamma"))
})

test_that("print.HRF reports its public summary and returns the object", {
  expect_output(
    returned <- print(HRF_SPMG3),
    "Basis functions: 3"
  )
  expect_identical(returned, HRF_SPMG3)
})

test_that("single HRF plots leave room for the peak label", {
  hrf <- gen_hrf(hrf_gamma, shape = 6, rate = 1)
  time <- seq(0, 25, by = 0.1)
  .with_pdf_device({
    plotted <- plot(hrf, time = time)
    expect_equal(plotted$response, hrf(time))
    expect_gt(graphics::par("usr")[4] - max(plotted$response),
              2 * graphics::strheight("Peak: 5.0 s", cex = 0.85))
    expect_no_error(plot(hrf, time = time, main = "Custom title", ylim = c(0, 0.3)))
  })
})

test_that("comparison plots can be customized without drawing", {
  skip_if_not_installed("ggplot2")
  original_device <- grDevices::dev.cur()
  time <- seq(0, 20, by = 0.5)
  h <- withVisible(plot_hrfs(HRF_SPMG1, HRF_GAMMA, time = time, draw = FALSE))
  r <- withVisible(plot_regressors(regressor(c(2, 8), HRF_SPMG1),
                                   grid = time, draw = FALSE))
  expect_identical(grDevices::dev.cur(), original_device)
  expect_false(h$visible)
  expect_false(r$visible)
  expect_equal(h$value$response, c(evaluate(HRF_SPMG1, time), evaluate(HRF_GAMMA, time)))
  expect_equal(r$value$response, evaluate(regressor(c(2, 8), HRF_SPMG1), time))
  expect_true(inherits(attr(h$value, "plot"), "ggplot"))
  expect_true(inherits(attr(r$value, "plot"), "ggplot"))
  expect_false(isTRUE(attr(attr(h$value, "plot")$theme, "complete")))
})

# capabilities("cairo") can be TRUE while the Cairo device cannot load (CRAN's
# macOS R links it against XQuartz). Find a PNG device type that actually
# writes a file without warnings.
.working_png_type <- function() {
  for (type in c("cairo", "quartz")) {
    path <- tempfile(fileext = ".png")
    ok <- tryCatch({
      grDevices::png(path, width = 50, height = 50, type = type)
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      grDevices::dev.off()
      file.exists(path) && file.info(path)$size > 0
    }, warning = function(w) FALSE, error = function(e) FALSE)
    if (!ok && grDevices::dev.cur() > 1 &&
        identical(names(grDevices::dev.cur()), "png")) {
      grDevices::dev.off()
    }
    unlink(path)
    if (isTRUE(ok)) return(type)
  }
  NULL
}

test_that("onset transparency is respected by event and feature plots", {
  png_type <- .working_png_type()
  skip_if(is.null(png_type), "no working PNG device")
  render <- function(x, show, alpha) {
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path), add = TRUE)
    grDevices::png(path, width = 600, height = 400, type = png_type)
    tryCatch(plot(x, grid = seq(0, 16, by = 0.1), show_onsets = show,
                  onset_alpha = alpha), finally = grDevices::dev.off())
    readBin(path, "raw", n = file.info(path)$size)
  }
  for (x in list(regressor(c(2, 8), HRF_SPMG1),
                 feature_regressor(c(1, 2, 1), dt = 2))) {
    hidden <- render(x, FALSE, 0)
    transparent <- render(x, TRUE, 0)
    opaque <- render(x, TRUE, 1)
    expect_identical(transparent, hidden)
    expect_false(identical(transparent, opaque))
  }
})
