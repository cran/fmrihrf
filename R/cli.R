#' Command line entrypoint for fmrihrf
#'
#' @param args Character vector of command line arguments.
#' @return Integer exit status. `0` indicates success, `1` indicates a domain
#'   failure, and `2` indicates usage or runtime errors.
#' @export
fmrihrf_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  tryCatch(
    .fmrihrf_cli_main(args),
    fmrihrf_cli_domain_error = function(e) {
      message(conditionMessage(e))
      1L
    },
    fmrihrf_cli_usage_error = function(e) {
      message(conditionMessage(e))
      2L
    },
    error = function(e) {
      message("fmrihrf: ", conditionMessage(e))
      2L
    }
  )
}

#' Install fmrihrf command line wrappers
#'
#' @param dest_dir Directory where wrapper commands should be copied.
#' @param overwrite Logical; overwrite an existing command only when `TRUE`.
#' @param commands Optional character vector of command names to install.
#' @return Invisibly, a named character vector of installed command paths.
#' @export
install_cli <- function(dest_dir = "~/.local/bin",
                        overwrite = FALSE,
                        commands = NULL) {
  command_map <- c(fmrihrf = "fmrihrf")
  if (is.null(commands)) {
    commands <- names(command_map)
  }
  unknown <- setdiff(commands, names(command_map))
  if (length(unknown) > 0) {
    stop("Unknown command(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }

  exec_dir <- .fmrihrf_exec_dir()
  if (!dir.exists(exec_dir)) {
    stop("Cannot find fmrihrf exec directory.", call. = FALSE)
  }

  dest_dir <- path.expand(dest_dir)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dest_dir)) {
    stop("Could not create destination directory: ", dest_dir, call. = FALSE)
  }

  installed <- stats::setNames(character(length(commands)), commands)
  for (cmd in commands) {
    src <- file.path(exec_dir, command_map[[cmd]])
    if (!file.exists(src)) {
      stop("Cannot find command wrapper: ", src, call. = FALSE)
    }
    dest <- file.path(dest_dir, cmd)
    if (file.exists(dest) && !isTRUE(overwrite)) {
      stop("Refusing to overwrite existing command: ", dest,
           "\nUse overwrite = TRUE to replace it.", call. = FALSE)
    }
    ok <- file.copy(src, dest, overwrite = TRUE)
    if (!ok) {
      stop("Could not copy command wrapper to: ", dest, call. = FALSE)
    }
    Sys.chmod(dest, mode = "0755")
    installed[[cmd]] <- dest
  }

  path_dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  if (!normalizePath(dest_dir, mustWork = FALSE) %in%
      normalizePath(path_dirs, mustWork = FALSE)) {
    message("Add ", dest_dir, " to PATH to run installed commands directly.")
  }

  invisible(installed)
}

.fmrihrf_cli_main <- function(args) {
  if (length(args) == 0 || args[[1]] %in% c("-h", "--help", "help")) {
    if (length(args) >= 2 && args[[1]] == "help") {
      return(.fmrihrf_cli_help(args[[2]]))
    }
    .fmrihrf_cli_help()
    return(0L)
  }

  command <- args[[1]]
  command_args <- args[-1]
  if ("--help" %in% command_args || "-h" %in% command_args) {
    return(.fmrihrf_cli_help(command))
  }

  switch(command,
    list = .fmrihrf_cli_list(command_args),
    eval = .fmrihrf_cli_eval(command_args),
    regressor = .fmrihrf_cli_regressor(command_args),
    design = .fmrihrf_cli_design(command_args),
    .cli_usage_error("Unknown command: ", command,
                     "\nRun `fmrihrf --help` for available commands.")
  )
}

.fmrihrf_cli_help <- function(command = NULL) {
  text <- switch(command %||% "main",
    main = c(
      "Usage: fmrihrf <command> [options]",
      "",
      "Commands:",
      "  list        List available HRFs and basis generators",
      "  eval        Evaluate an HRF or basis over a time grid",
      "  regressor   Build and evaluate one event regressor",
      "  design      Build a condition design matrix from an events table",
      "",
      "Run `fmrihrf help <command>` for command-specific options."
    ),
    list = c(
      "Usage: fmrihrf list [--details] [--json] [--output FILE]",
      "",
      "Options:",
      "  --details       Include descriptions from the HRF registry",
      "  --json          Write JSON instead of CSV",
      "  --output FILE   Write output to FILE instead of stdout"
    ),
    eval = c(
      "Usage: fmrihrf eval [options]",
      "",
      "Core options:",
      "  --hrf NAME          HRF name from `fmrihrf list` [spmg1]",
      "  --from SEC          Grid start in seconds [0]",
      "  --to SEC            Grid end in seconds [32]",
      "  --by SEC            Grid step in seconds [0.5]",
      "  --times LIST        Comma-separated grid values; overrides from/to/by",
      "  --nbasis N          Basis count for generated HRFs [5]",
      "  --span SEC          HRF span in seconds [24]",
      "  --lag SEC           Shift HRF by this many seconds [0]",
      "  --width SEC         Block HRF width in seconds [0]",
      "  --amplitude VALUE   Event amplitude for HRF evaluation [1]",
      "  --duration SEC      Event duration for HRF evaluation [0]",
      "  --precision SEC     Internal block precision [0.2]",
      "  --normalize         Normalize columns to unit peak",
      "  --no-summate        Average block response instead of summating",
      "  --json              Write JSON instead of CSV",
      "  --output FILE       Write output to FILE instead of stdout"
    ),
    regressor = c(
      "Usage: fmrihrf regressor [options]",
      "",
      "Options:",
      "  --onsets LIST       Comma-separated onset times in seconds",
      "  --events FILE       CSV/TSV with onset[,duration,amplitude] columns",
      "  --hrf NAME          HRF name from `fmrihrf list` [spmg1]",
      "  --blocklens LIST    Scans per block for acquisition grid",
      "  --tr LIST           Repetition time(s) in seconds",
      "  --from SEC          Grid start when no sampling frame is supplied [0]",
      "  --to SEC            Grid end when no sampling frame is supplied",
      "  --by SEC            Grid step when no sampling frame is supplied [1]",
      "  --duration LIST     Event duration(s) when --events is not used [0]",
      "  --amplitude LIST    Event amplitude(s) when --events is not used [1]",
      "  --nbasis N          Basis count for generated HRFs [5]",
      "  --span SEC          HRF/regressor span in seconds [24]",
      "  --precision SEC     Internal convolution precision [0.33]",
      "  --method METHOD     conv, fft, Rconv, or loop [conv]",
      "  --normalize         Normalize output columns to unit peak",
      "  --json              Write JSON instead of CSV",
      "  --output FILE       Write output to FILE instead of stdout"
    ),
    design = c(
      "Usage: fmrihrf design --events FILE --blocklens LIST --tr LIST [options]",
      "",
      "The events file must include onset, condition, and block columns.",
      "Optional duration and amplitude columns are used when present.",
      "",
      "Options:",
      "  --events FILE       CSV/TSV events table",
      "  --condition COL     Condition column name [condition]",
      "  --onset COL         Onset column name [onset]",
      "  --block COL         Block column name [block]",
      "  --duration COL      Duration column name [duration]",
      "  --amplitude COL     Amplitude column name [amplitude]",
      "  --blocklens LIST    Scans per block",
      "  --tr LIST           Repetition time(s) in seconds",
      "  --start-time LIST   First acquisition offset(s); default TR / 2",
      "  --hrf NAME          HRF name from `fmrihrf list` [spmg1]",
      "  --nbasis N          Basis count for generated HRFs [5]",
      "  --span SEC          HRF/regressor span in seconds [24]",
      "  --precision SEC     Internal convolution precision [0.33]",
      "  --method METHOD     conv, fft, Rconv, or loop [conv]",
      "  --sparse            Return sparse design internally before writing dense output",
      "  --json              Write JSON instead of CSV",
      "  --output FILE       Write output to FILE instead of stdout"
    ),
    .cli_usage_error("Unknown help topic: ", command)
  )
  cat(paste(text, collapse = "\n"), "\n", sep = "")
  0L
}

.fmrihrf_cli_list <- function(args) {
  opts <- .parse_cli_args(args, flags = c("details", "json"),
                         values = c("output"))
  info <- list_available_hrfs(details = isTRUE(opts$details))
  .write_cli_table(info, json = isTRUE(opts$json), output = opts$output)
  0L
}

.fmrihrf_cli_eval <- function(args) {
  opts <- .parse_cli_args(
    args,
    flags = c("normalize", "json", "summate"),
    false_flags = c("summate"),
    values = c("hrf", "from", "to", "by", "times", "nbasis", "span", "lag",
               "width", "amplitude", "duration", "precision", "output")
  )
  opts <- .defaults(opts, list(
    hrf = "spmg1", from = "0", to = "32", by = "0.5",
    nbasis = "5", span = "24", lag = "0", width = "0",
    amplitude = "1", duration = "0", precision = "0.2",
    summate = TRUE
  ))

  grid <- .grid_from_options(opts)
  hrf <- getHRF(
    opts$hrf,
    nbasis = .as_scalar_integer(opts$nbasis, "nbasis"),
    span = .as_scalar_numeric(opts$span, "span"),
    lag = .as_scalar_numeric(opts$lag, "lag"),
    width = .as_scalar_numeric(opts$width, "width"),
    summate = isTRUE(opts$summate),
    normalize = isTRUE(opts$normalize)
  )
  values <- evaluate(
    hrf, grid,
    amplitude = .as_scalar_numeric(opts$amplitude, "amplitude"),
    duration = .as_scalar_numeric(opts$duration, "duration"),
    precision = .as_scalar_numeric(opts$precision, "precision"),
    summate = isTRUE(opts$summate),
    normalize = isTRUE(opts$normalize)
  )

  out <- .values_to_frame(grid, values, time_name = "time")
  .write_cli_table(out, json = isTRUE(opts$json), output = opts$output)
  0L
}

.fmrihrf_cli_regressor <- function(args) {
  opts <- .parse_cli_args(
    args,
    flags = c("normalize", "json"),
    values = c("onsets", "events", "hrf", "blocklens", "tr", "start-time",
               "from", "to", "by", "duration", "amplitude", "nbasis", "span",
               "precision", "method", "output")
  )
  opts <- .defaults(opts, list(
    hrf = "spmg1", by = "1", duration = "0", amplitude = "1",
    nbasis = "5", span = "24", precision = "0.33", method = "conv"
  ))

  events <- .events_for_regressor(opts)
  hrf <- getHRF(opts$hrf,
                nbasis = .as_scalar_integer(opts$nbasis, "nbasis"),
                span = .as_scalar_numeric(opts$span, "span"))
  reg <- regressor(
    onsets = events$onset,
    hrf = hrf,
    duration = events$duration,
    amplitude = events$amplitude,
    span = .as_scalar_numeric(opts$span, "span")
  )
  grid <- .regressor_grid(opts, events$onset)
  values <- evaluate(
    reg, grid = grid,
    precision = .as_scalar_numeric(opts$precision, "precision"),
    method = .match_choice(opts$method, c("conv", "loop", "fft", "Rconv"), "method"),
    normalize = isTRUE(opts$normalize)
  )

  out <- .values_to_frame(grid, values, time_name = "time")
  .write_cli_table(out, json = isTRUE(opts$json), output = opts$output)
  0L
}

.fmrihrf_cli_design <- function(args) {
  opts <- .parse_cli_args(
    args,
    flags = c("sparse", "json"),
    values = c("events", "condition", "onset", "block", "duration", "amplitude",
               "blocklens", "tr", "start-time", "hrf", "nbasis", "span",
               "precision", "method", "output")
  )
  opts <- .defaults(opts, list(
    condition = "condition", onset = "onset", block = "block",
    duration = "duration", amplitude = "amplitude", hrf = "spmg1",
    nbasis = "5", span = "24", precision = "0.33", method = "conv"
  ))
  .require_options(opts, c("events", "blocklens", "tr"))

  events <- .read_events_table(opts$events)
  .require_columns(events, c(opts$onset, opts$condition, opts$block))
  duration <- if (opts$duration %in% names(events)) events[[opts$duration]] else 0
  amplitude <- if (opts$amplitude %in% names(events)) events[[opts$amplitude]] else 1
  sframe <- .sampling_frame_from_options(opts)
  hrf <- getHRF(opts$hrf,
                nbasis = .as_scalar_integer(opts$nbasis, "nbasis"),
                span = .as_scalar_numeric(opts$span, "span"))

  design <- regressor_design(
    onsets = events[[opts$onset]],
    fac = events[[opts$condition]],
    block = events[[opts$block]],
    sframe = sframe,
    hrf = hrf,
    duration = duration,
    amplitude = amplitude,
    span = .as_scalar_numeric(opts$span, "span"),
    precision = .as_scalar_numeric(opts$precision, "precision"),
    method = .match_choice(opts$method, c("conv", "loop", "fft", "Rconv"), "method"),
    sparse = isTRUE(opts$sparse)
  )

  grid <- samples(sframe, global = TRUE)
  design <- as.matrix(design)
  design_names <- .design_column_names(events[[opts$condition]], nbasis(hrf))
  if (length(design_names) == ncol(design)) {
    colnames(design) <- design_names
  }
  out <- .values_to_frame(grid, design, time_name = "time")
  .write_cli_table(out, json = isTRUE(opts$json), output = opts$output)
  0L
}

.parse_cli_args <- function(args, flags = character(), values = character(),
                            false_flags = character()) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      .cli_usage_error("Unexpected positional argument: ", token)
    }
    token <- substring(token, 3L)
    value <- NULL
    if (grepl("=", token, fixed = TRUE)) {
      parts <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else {
      key <- token
    }

    if (startsWith(key, "no-")) {
      flag <- substring(key, 4L)
      if (!flag %in% false_flags) {
        .cli_usage_error("Unknown option: --", key)
      }
      if (!is.null(value)) {
        .cli_usage_error("Boolean option does not take a value: --", key)
      }
      out[[flag]] <- FALSE
      i <- i + 1L
      next
    }

    if (key %in% flags) {
      if (!is.null(value)) {
        .cli_usage_error("Boolean option does not take a value: --", key)
      }
      out[[key]] <- TRUE
      i <- i + 1L
      next
    }

    if (key %in% values) {
      if (is.null(value)) {
        if (i == length(args)) {
          .cli_usage_error("Missing value for option: --", key)
        }
        value <- args[[i + 1L]]
        if (startsWith(value, "--")) {
          .cli_usage_error("Missing value for option: --", key)
        }
        i <- i + 1L
      }
      out[[key]] <- value
      i <- i + 1L
      next
    }

    .cli_usage_error("Unknown option: --", key)
  }
  out
}

.defaults <- function(opts, defaults) {
  for (nm in names(defaults)) {
    if (is.null(opts[[nm]])) {
      opts[[nm]] <- defaults[[nm]]
    }
  }
  opts
}

.require_options <- function(opts, names) {
  missing <- names[vapply(names, function(nm) is.null(opts[[nm]]), logical(1))]
  if (length(missing) > 0) {
    .cli_usage_error("Missing required option(s): ",
                     paste(paste0("--", missing), collapse = ", "))
  }
}

.require_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    .cli_domain_error("Missing required column(s): ", paste(missing, collapse = ", "))
  }
}

.grid_from_options <- function(opts) {
  if (!is.null(opts$times)) {
    values <- .parse_numeric_list(opts$times, "times")
    if (length(values) == 0) {
      .cli_usage_error("--times must contain at least one number")
    }
    return(values)
  }
  from <- .as_scalar_numeric(opts$from, "from")
  to <- .as_scalar_numeric(opts$to, "to")
  by <- .as_scalar_numeric(opts$by, "by")
  if (by <= 0) {
    .cli_usage_error("--by must be positive")
  }
  if (to < from) {
    .cli_usage_error("--to must be greater than or equal to --from")
  }
  seq(from, to, by = by)
}

.regressor_grid <- function(opts, onsets) {
  if (!is.null(opts$blocklens) || !is.null(opts$tr)) {
    .require_options(opts, c("blocklens", "tr"))
    return(samples(.sampling_frame_from_options(opts), global = TRUE))
  }
  if (is.null(opts$from)) {
    opts$from <- "0"
  }
  if (is.null(opts$to)) {
    max_onset <- if (length(onsets) > 0) max(onsets) else 0
    opts$to <- as.character(max_onset + .as_scalar_numeric(opts$span, "span"))
  }
  .grid_from_options(opts)
}

.sampling_frame_from_options <- function(opts) {
  blocklens <- .parse_numeric_list(opts$blocklens, "blocklens")
  tr <- .parse_numeric_list(opts$tr, "tr")
  start_time <- if (is.null(opts[["start-time"]])) {
    tr / 2
  } else {
    .parse_numeric_list(opts[["start-time"]], "start-time")
  }
  sampling_frame(
    blocklens = blocklens,
    TR = tr,
    start_time = start_time,
    precision = .as_scalar_numeric(opts$precision, "precision")
  )
}

.events_for_regressor <- function(opts) {
  if (!is.null(opts$events)) {
    events <- .read_events_table(opts$events)
    .require_columns(events, "onset")
    duration <- if ("duration" %in% names(events)) events$duration else 0
    amplitude <- if ("amplitude" %in% names(events)) events$amplitude else 1
    return(data.frame(
      onset = as.numeric(events$onset),
      duration = as.numeric(duration),
      amplitude = as.numeric(amplitude)
    ))
  }
  if (is.null(opts$onsets)) {
    .cli_usage_error("Provide --onsets or --events")
  }
  onset <- .parse_numeric_list(opts$onsets, "onsets")
  duration <- recycle_or_error(.parse_numeric_list(opts$duration, "duration"),
                               length(onset), "duration")
  amplitude <- recycle_or_error(.parse_numeric_list(opts$amplitude, "amplitude"),
                                length(onset), "amplitude")
  data.frame(
    onset = onset,
    duration = duration,
    amplitude = amplitude
  )
}

.read_events_table <- function(path) {
  if (!file.exists(path)) {
    .cli_domain_error("File does not exist: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "tab")) "\t" else ","
  utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "", quote = "\"")
}

.parse_numeric_list <- function(x, name) {
  if (is.numeric(x)) {
    return(x)
  }
  if (length(x) != 1 || is.na(x)) {
    .cli_usage_error("--", name, " must be a scalar or comma-separated list")
  }
  pieces <- strsplit(x, ",", fixed = TRUE)[[1]]
  pieces <- trimws(pieces)
  if (any(!nzchar(pieces))) {
    .cli_usage_error("--", name, " contains an empty value")
  }
  values <- suppressWarnings(as.numeric(pieces))
  if (anyNA(values)) {
    .cli_usage_error("--", name, " must contain only numeric values")
  }
  values
}

.as_scalar_numeric <- function(x, name) {
  values <- .parse_numeric_list(x, name)
  if (length(values) != 1) {
    .cli_usage_error("--", name, " must be a single number")
  }
  values[[1]]
}

.as_scalar_integer <- function(x, name) {
  value <- .as_scalar_numeric(x, name)
  if (!is.finite(value) || value %% 1 != 0) {
    .cli_usage_error("--", name, " must be a whole number")
  }
  as.integer(value)
}

.match_choice <- function(value, choices, name) {
  if (!value %in% choices) {
    .cli_usage_error("--", name, " must be one of: ",
                     paste(choices, collapse = ", "))
  }
  value
}

.values_to_frame <- function(grid, values, time_name = "time") {
  if (inherits(values, "Matrix")) {
    values <- as.matrix(values)
  }
  if (is.null(dim(values))) {
    out <- data.frame(grid, value = as.numeric(values), check.names = FALSE)
  } else {
    values <- as.matrix(values)
    colnames(values) <- colnames(values) %||% paste0("basis", seq_len(ncol(values)))
    empty <- !nzchar(colnames(values))
    colnames(values)[empty] <- paste0("basis", which(empty))
    out <- data.frame(grid, values, check.names = FALSE)
  }
  names(out)[[1]] <- time_name
  out
}

.design_column_names <- function(condition, nbasis) {
  levels <- levels(as.factor(condition))
  if (identical(as.integer(nbasis), 1L)) {
    return(levels)
  }
  unlist(lapply(levels, function(level) {
    paste0(level, "_basis", seq_len(nbasis))
  }), use.names = FALSE)
}

.write_cli_table <- function(data, json = FALSE, output = NULL) {
  if (isTRUE(json)) {
    text <- jsonlite::toJSON(data, dataframe = "rows", auto_unbox = TRUE,
                             pretty = TRUE, na = "null")
  } else {
    buffer <- character()
    con <- textConnection("buffer", "w", local = TRUE)
    on.exit(close(con), add = TRUE)
    utils::write.csv(data, con, row.names = FALSE, na = "")
    text <- paste(buffer, collapse = "\n")
  }

  if (is.null(output)) {
    cat(text, "\n", sep = "")
  } else {
    writeLines(text, output, useBytes = TRUE)
  }
}

.cli_usage_error <- function(...) {
  msg <- paste0(...)
  stop(structure(list(message = msg, call = NULL),
                 class = c("fmrihrf_cli_usage_error", "error", "condition")))
}

.cli_domain_error <- function(...) {
  msg <- paste0(...)
  stop(structure(list(message = msg, call = NULL),
                 class = c("fmrihrf_cli_domain_error", "error", "condition")))
}

.fmrihrf_exec_dir <- function() {
  exec_dir <- system.file("exec", package = "fmrihrf")
  if (nzchar(exec_dir) && dir.exists(exec_dir)) {
    return(exec_dir)
  }

  source_exec <- file.path(getwd(), "exec")
  if (dir.exists(source_exec)) {
    return(source_exec)
  }

  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkg_path <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
    if (!is.null(pkg_path)) {
      pkg_exec <- file.path(pkg_path, "exec")
      if (dir.exists(pkg_exec)) {
        return(pkg_exec)
      }
    }
  }

  ""
}
