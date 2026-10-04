test_that("fmrihrf_cli help paths return success", {
  out <- capture.output(status <- fmrihrf_cli("--help"))
  expect_equal(status, 0L)
  expect_true(any(grepl("Usage: fmrihrf", out, fixed = TRUE)))

  out <- capture.output(status <- fmrihrf_cli(c("help", "eval")))
  expect_equal(status, 0L)
  expect_true(any(grepl("Usage: fmrihrf eval", out, fixed = TRUE)))

  out <- capture.output(status <- fmrihrf_cli(c("eval", "--help")))
  expect_equal(status, 0L)
  expect_true(any(grepl("Usage: fmrihrf eval", out, fixed = TRUE)))
})

test_that("parser rejects unknown options and missing values", {
  expect_equal(suppressMessages(fmrihrf_cli(c("list", "--bogus"))), 2L)
  expect_equal(suppressMessages(fmrihrf_cli(c("eval", "--hrf"))), 2L)
})

test_that("parser supports equals values and negative numeric values", {
  out <- capture.output(
    status <- fmrihrf_cli(c("eval", "--hrf=spmg1", "--times", "-0.5,0,0.5"))
  )
  expect_equal(status, 0L)
  expect_equal(out[1], "\"time\",\"value\"")
  expect_match(out[2], "^-0.5,")
})

test_that("list command emits parseable JSON", {
  out <- capture.output(status <- fmrihrf_cli(c("list", "--json")))
  expect_equal(status, 0L)
  parsed <- jsonlite::fromJSON(paste(out, collapse = "\n"))
  expect_true(all(c("name", "type", "nbasis_default") %in% names(parsed)))
  expect_true("spmg1" %in% parsed$name)
})

test_that("eval command emits expected columns", {
  path <- tempfile(fileext = ".csv")
  status <- fmrihrf_cli(c(
    "eval", "--hrf", "spmg3", "--from", "0", "--to", "1", "--by", "1",
    "--output", path
  ))
  expect_equal(status, 0L)
  out <- utils::read.csv(path, check.names = FALSE)
  expect_equal(names(out), c("time", "basis1", "basis2", "basis3"))
  expect_equal(out$time, c(0, 1))
})

test_that("regressor command evaluates onsets over a sampling frame", {
  path <- tempfile(fileext = ".csv")
  status <- fmrihrf_cli(c(
    "regressor", "--onsets", "0,4", "--blocklens", "4", "--tr", "1",
    "--hrf", "spmg1", "--output", path
  ))
  expect_equal(status, 0L)
  out <- utils::read.csv(path, check.names = FALSE)
  expect_equal(names(out), c("time", "value"))
  expect_equal(nrow(out), 4L)
})

test_that("design command creates condition columns from an events table", {
  events <- tempfile(fileext = ".csv")
  writeLines(c(
    "onset,condition,block,duration,amplitude",
    "0,A,1,0,1",
    "4,B,1,0,1"
  ), events)
  path <- tempfile(fileext = ".csv")
  status <- fmrihrf_cli(c(
    "design", "--events", events, "--blocklens", "6", "--tr", "1",
    "--output", path
  ))
  expect_equal(status, 0L)
  out <- utils::read.csv(path, check.names = FALSE)
  expect_equal(names(out), c("time", "A", "B"))
  expect_equal(nrow(out), 6L)
})

test_that("install_cli copies the wrapper and refuses accidental overwrite", {
  dest_dir <- tempfile()
  installed <- install_cli(dest_dir)
  expect_true(file.exists(installed[["fmrihrf"]]))
  if (.Platform$OS.type != "windows") {
    expect_equal(unname(file.access(installed[["fmrihrf"]], mode = 1)), 0L)
  }
  expect_error(install_cli(dest_dir), "Refusing to overwrite")
})

test_that("exec wrapper has a valid help path", {
  wrapper <- system.file("exec", "fmrihrf", package = "fmrihrf")
  if (!nzchar(wrapper)) {
    wrapper <- file.path(getwd(), "exec", "fmrihrf")
  }
  skip_if_not(file.exists(wrapper))
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  result <- system2(rscript,
                    c(shQuote(wrapper), "--help"), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (is.null(status)) {
    status <- 0L
  }
  expect_equal(status, 0L)
  expect_true(any(grepl("Usage: fmrihrf", result, fixed = TRUE)))
})
