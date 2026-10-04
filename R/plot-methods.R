#' Plot an HRF Object
#'
#' Draws an HRF with base graphics. Single-basis HRFs show the response curve
#' with its peak annotated. Multi-basis HRFs (e.g., [HRF_SPMG3]) show every
#' basis function, coloured along the ordered [hrf_palette()].
#'
#' @md
#' @param x An HRF object
#' @param time Numeric vector of time points. If NULL (default), uses
#'   seq(0, span, by = 0.1) where span is the HRF's span attribute.
#' @param normalize Logical; if TRUE, normalize responses to peak at 1.
#'   Default is FALSE.
#' @param show_peak Logical; if TRUE (default for single-basis HRFs), annotate
#'   the peak time on the plot.
#' @param ... Additional arguments passed to [graphics::plot()], such as
#'   `main` or `ylim`.
#' @return Invisibly returns a data frame with the time and response values
#'   (useful for further customization).
#' @seealso [plot_hrfs()] for ggplot2 comparisons of several HRFs.
#' @examples
#' # Plot single-basis HRF
#' plot(HRF_SPMG1)
#'
#' # Plot multi-basis HRF
#' plot(HRF_SPMG3)
#'
#' # Plot with normalization
#' plot(HRF_GAMMA, normalize = TRUE)
#'
#' # Custom time range
#' plot(HRF_SPMG1, time = seq(0, 30, by = 0.5))
#' @method plot HRF
#' @export
plot.HRF <- function(x, time = NULL, normalize = FALSE, show_peak = TRUE, ...) {
  span <- attr(x, "span") %||% 24
  if (is.null(time)) {
    time <- seq(0, span, by = 0.1)
  }

  y <- evaluate(x, time, normalize = normalize)
  hrf_name <- attr(x, "name") %||% "HRF"
  args <- list(...)
  main <- args$main %||% hrf_name
  args$main <- NULL
  ylab <- if (normalize) "Response / peak" else "Response"

  if (is.matrix(y)) {
    nb <- ncol(y)
    cols <- .basis_palette(nb)
    do.call(.base_series_panel, c(list(time = time, Y = y, colours = cols,
                                       labels = paste0("B", seq_len(nb)),
                                       main = main, ylab = ylab), args))
    df <- data.frame(time = time, y)
    colnames(df)[-1] <- paste0("basis_", seq_len(nb))
  } else {
    if (show_peak && length(y) > 0 && is.null(args$ylim)) {
      limits <- range(y, 0, finite = TRUE)
      limits[2] <- limits[2] + max(diff(limits), abs(limits[2]), 1e-6) * 0.18
      args$ylim <- limits
    }
    do.call(.base_series_panel, c(list(time = time, Y = y,
                                       colours = hrf_palette(1),
                                       main = main, ylab = ylab), args))
    if (show_peak && length(y) > 0) {
      peak_idx <- which.max(y)
      graphics::points(time[peak_idx], y[peak_idx], pch = 19, cex = 0.9,
                       col = hrf_palette(1))
      graphics::text(time[peak_idx], y[peak_idx],
                     sprintf("Peak: %.1f s", time[peak_idx]),
                     pos = 3, offset = 0.5, cex = 0.85, col = .hrf_neutral)
    }
    df <- data.frame(time = time, response = y)
  }

  invisible(df)
}

#' Compare Multiple HRF Functions
#'
#' Plots one or more HRF objects on shared axes. Multi-basis HRFs are
#' expanded into one curve per basis function, so `plot_hrfs(HRF_SPMG3)`
#' shows the canonical response and both derivatives. Uses ggplot2 when
#' available, otherwise base graphics. Colours come from [hrf_palette()].
#'
#' Inside a knitr document the result is returned visibly and printed by
#' knitr, which lets document themes (for example dark-mode figure twins)
#' handle the ggplot. At the console the plot is drawn immediately and the
#' data are returned invisibly.
#'
#' @md
#' @param ... HRF objects to compare. Can be passed as individual arguments
#'   or as a named list.
#' @param time Numeric vector of time points. If NULL (default), uses
#'   seq(0, max_span, by = 0.1) where max_span is the maximum span across
#'   all HRFs.
#' @param normalize Logical; if TRUE, normalize all HRFs to peak at 1.
#'   Useful for comparing shapes regardless of amplitude. Default is FALSE.
#' @param labels Character vector of labels, either one per HRF or one per
#'   plotted curve (after basis expansion). If NULL (default), uses the 'name'
#'   attribute of each HRF; basis functions of a single basis set are labelled
#'   `B1`, `B2`, and so on.
#' @param title Character string for the plot title. If NULL (default), uses
#'   the HRF name for a single HRF and "HRF comparison" otherwise.
#' @param subtitle Character string for the plot subtitle. If NULL (default),
#'   no subtitle is shown.
#' @param use_ggplot Logical; if TRUE and ggplot2 is available, use ggplot2
#'   for plotting. If FALSE, use base R graphics. Default is TRUE.
#' @param draw Logical; draw the plot (default TRUE). With `use_ggplot = TRUE`,
#'   set FALSE to customize the returned data frame's `"plot"` attribute.
#' @param basis Either `"all"` (default) to plot every basis function of a
#'   multi-basis HRF, or `"first"` to plot only its first column.
#' @param layout Either `"overlay"` (default) to draw all curves in one panel
#'   or `"stack"` to give each curve its own panel on a shared time axis.
#' @param scales For `layout = "stack"`: `"free_y"` (default) or `"fixed"`
#'   (one shared y range).
#' @param palette Colour palette: `"auto"` (default) uses the ordered palette
#'   for a single basis set with more than three functions and the
#'   categorical palette otherwise; `"categorical"` or `"ordered"` force one
#'   (see [hrf_palette()]). Use `"ordered"` when the HRFs differ along one
#'   parameter (lag, width).
#' @param reference Optional HRF drawn as a dashed grey reference curve (for
#'   example the canonical HRF), outside the colour legend. It is evaluated on
#'   the same time grid and normalized like the other HRFs.
#' @param reference_label Text naming the reference curve, shown as a caption.
#'   Defaults to the reference HRF's name.
#' @return A data frame in long format with columns 'time', 'HRF', and
#'   'response'. With ggplot2, the plot is stored in the `"plot"` attribute.
#'   The value is invisible except when drawing inside knitr, where it prints
#'   like a ggplot object: it is drawn when it is the visible value of a chunk
#'   (not when assigned), and subsetting or modifying it, or calling
#'   `as.data.frame()`, returns a plain data frame. Call `as.data.frame()`
#'   before combining results with `rbind()` or dplyr verbs.
#' @examples
#' # Compare canonical HRFs
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, HRF_GAUSSIAN)
#'
#' # A basis set: one curve per basis function
#' plot_hrfs(HRF_SPMG3,
#'           labels = c("Canonical", "Temporal derivative", "Dispersion derivative"))
#'
#' # HRFs ordered by a parameter use the ordered palette
#' plot_hrfs(block_hrf(HRF_SPMG1, width = 1), block_hrf(HRF_SPMG1, width = 3),
#'           block_hrf(HRF_SPMG1, width = 5),
#'           labels = c("1 s", "3 s", "5 s"), palette = "ordered",
#'           title = "Effect of event duration")
#'
#' # Normalize for shape comparison
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, HRF_GAUSSIAN, normalize = TRUE,
#'           subtitle = "All HRFs normalized to peak at 1")
#'
#' # Use base R graphics instead of ggplot2
#' plot_hrfs(HRF_SPMG1, HRF_GAMMA, use_ggplot = FALSE)
#' @export
plot_hrfs <- function(..., time = NULL, normalize = FALSE, labels = NULL,
                      title = NULL, subtitle = NULL, use_ggplot = TRUE,
                      draw = TRUE, basis = c("all", "first"),
                      layout = c("overlay", "stack"),
                      palette = c("auto", "categorical", "ordered"),
                      reference = NULL, reference_label = NULL,
                      scales = c("free_y", "fixed")) {
  scales <- match.arg(scales)
  basis <- match.arg(basis)
  layout <- match.arg(layout)
  palette <- match.arg(palette)
  if (!is.null(reference) && !inherits(reference, "HRF")) {
    stop("'reference' must be an HRF object")
  }
  hrfs <- list(...)

  # Handle case where a single list is passed
  if (length(hrfs) == 1 && is.list(hrfs[[1]]) && !inherits(hrfs[[1]], "HRF")) {
    hrfs <- hrfs[[1]]
  }

  n_hrfs <- length(hrfs)
  if (n_hrfs == 0) {
    stop("At least one HRF object must be provided")
  }
  for (i in seq_along(hrfs)) {
    if (!inherits(hrfs[[i]], "HRF")) {
      stop("All arguments must be HRF objects. Argument ", i, " is not an HRF.")
    }
  }

  if (is.null(time)) {
    spans <- vapply(hrfs, function(h) as.numeric(attr(h, "span") %||% 24), numeric(1))
    time <- seq(0, max(spans), by = 0.1)
  }

  # Evaluate, expanding multi-basis HRFs into columns.
  mats <- lapply(hrfs, function(h) {
    y <- as.matrix(evaluate(h, time, normalize = normalize))
    if (basis == "first") y[, 1, drop = FALSE] else y
  })
  nb <- vapply(mats, ncol, integer(1))
  n_series <- sum(nb)
  single_set <- n_hrfs == 1 && nb[1] > 1

  hrf_names <- vapply(hrfs, function(h) as.character(attr(h, "name") %||% "HRF")[1],
                      character(1))
  if (is.null(labels)) {
    base <- hrf_names
    if (any(duplicated(base))) base <- paste0(base, "_", seq_along(base))
  } else if (length(labels) == n_series) {
    base <- NULL
  } else if (length(labels) == n_hrfs) {
    base <- labels
  } else {
    stop("Length of 'labels' must match number of HRFs")
  }
  if (is.null(base)) {
    series_labels <- as.character(labels)
  } else {
    series_labels <- unlist(lapply(seq_len(n_hrfs), function(i) {
      if (nb[i] == 1) return(base[i])
      if (single_set) paste0("B", seq_len(nb[i])) else paste0(base[i], " B", seq_len(nb[i]))
    }))
  }
  if (any(duplicated(series_labels))) {
    series_labels <- make.unique(series_labels, sep = " ")
  }

  resp <- do.call(cbind, mats)
  df <- data.frame(
    time = rep(time, n_series),
    HRF = factor(rep(series_labels, each = length(time)), levels = series_labels),
    response = as.vector(resp)
  )

  if (is.null(title)) {
    title <- if (n_hrfs == 1) hrf_names[1] else "HRF comparison"
  }
  pal <- .series_palette(n_series, palette, ordered_hint = single_set && nb[1] > 3)
  ylab <- if (normalize) "Response / peak" else "Response"
  ref_df <- NULL
  if (!is.null(reference)) {
    ref_y <- as.matrix(evaluate(reference, time, normalize = normalize))[, 1]
    ref_df <- data.frame(time = time, response = ref_y)
    reference_label <- reference_label %||% as.character(attr(reference, "name") %||% "reference")
  }

  if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    gdf <- df
    names(gdf)[2] <- "series"
    direct <- single_set && is.null(labels) && layout == "overlay"
    p <- .gg_series(gdf, pal$colours, title = title, subtitle = subtitle,
                    ylab = ylab, layout = layout, direct_labels = direct,
                    linetypes = pal$type == "categorical" && n_series > 4,
                    reference = ref_df, scales = scales)
    if (!is.null(ref_df)) {
      p <- p + ggplot2::labs(caption = paste("Dashed grey:", reference_label)) +
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0),
                       plot.caption.position = "plot")
    }
    return(.finish_gg(df, p, draw))
  }
  if (draw) {
    .base_series_panel(time, resp, pal$colours, labels = series_labels,
                       main = title, ylab = ylab)
    if (!is.null(ref_df)) {
      graphics::lines(ref_df$time, ref_df$response, col = .hrf_neutral, lty = 2)
    }
  }
  invisible(df)
}

#' Plot a Regressor Object
#'
#' Draws the predicted BOLD time course of a regressor with base graphics.
#' Event onsets are marked with ticks along the time axis. Regressors built
#' from a basis set are drawn one basis function per panel by default.
#'
#' @md
#' @param x A `Reg` object created by `regressor()`.
#' @param grid Numeric vector of time points for evaluation. If NULL (default),
#'   uses a grid from 0 to max(onsets) + span with step 0.25 s.
#' @param show_onsets Logical; if TRUE (default), mark event onsets with ticks
#'   on the time axis.
#' @param onset_color Colour for onset ticks. If NULL (default), a neutral grey.
#' @param onset_alpha Alpha transparency for onset ticks. Default is 0.5.
#' @param precision Numeric sampling precision for HRF evaluation. If NULL
#'   (default), the grid spacing capped at 0.33 s, so sharp HRF edges are
#'   drawn where they occur.
#' @param layout For multi-basis regressors, `"stack"` (default) draws one
#'   panel per basis function; `"overlay"` draws them in a single panel.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return Invisibly returns a data frame with the time and response values.
#' @seealso [plot_regressors()] for ggplot2 comparisons of several regressors.
#' @examples
#' # Create and plot a simple regressor
#' reg <- regressor(onsets = c(10, 30, 50), hrf = HRF_SPMG1)
#' plot(reg)
#'
#' # Plot with custom time grid
#' plot(reg, grid = seq(0, 80, by = 1))
#'
#' # Plot without onset markers
#' plot(reg, show_onsets = FALSE)
#'
#' # A basis-set regressor: one panel per basis function
#' plot(regressor(c(10, 40), HRF_SPMG3))
#' @method plot Reg
#' @export
plot.Reg <- function(x, grid = NULL, show_onsets = TRUE,
                     onset_color = NULL, onset_alpha = 0.5,
                     precision = NULL, layout = c("stack", "overlay"), ...) {
  layout <- match.arg(layout)
  if (is.null(grid)) {
    max_time <- max(x$onsets, na.rm = TRUE) + x$span
    grid <- seq(0, max_time, by = 0.25)
  }
  response <- evaluate(x, grid, precision = .plot_precision(precision, grid))

  hrf_is_list <- isTRUE(attr(x, "hrf_is_list"))
  hrf_name <- if (hrf_is_list) {
    "trial-varying HRF"
  } else {
    attr(x$hrf, "name") %||% "custom HRF"
  }
  n_events <- length(x$onsets)
  title <- sprintf("%s regressor, %d event%s", hrf_name, n_events,
                   if (n_events == 1) "" else "s")
  .plot_reg_base(response, grid, x$onsets, title, show_onsets, onset_color,
                 onset_alpha, layout, ...)
}

# Shared base-graphics body for plot.Reg() and plot.FeatureReg().
.plot_reg_base <- function(response, grid, onsets, title, show_onsets,
                           onset_color, onset_alpha, layout, ...) {
  onset_col <- grDevices::adjustcolor(onset_color %||% .hrf_neutral,
                                      alpha.f = onset_alpha)
  marks <- if (show_onsets) onsets else NULL

  if (is.matrix(response)) {
    nb <- ncol(response)
    cols <- .basis_palette(nb)
    if (layout == "stack" && nb > 1) {
      op <- graphics::par(mfrow = c(nb, 1), mar = c(1.6, 4.1, 1.4, 1.1),
                          oma = c(2.6, 0, 1.8, 0))
      on.exit(graphics::par(op), add = TRUE)
      for (j in seq_len(nb)) {
        .base_series_panel(grid, response[, j], cols[j], main = NULL,
                           xlab = "", ylab = paste0("B", j), onsets = marks,
                           onset_col = onset_col, legend = FALSE,
                           xaxt = if (j < nb) "n" else "s", ...)
      }
      graphics::mtext("Time (s)", side = 1, outer = TRUE, line = 1.2,
                      cex = graphics::par("cex.lab") %||% 1)
      graphics::mtext(title, side = 3, outer = TRUE, line = 0.4, adj = 0,
                      font = 2)
    } else {
      .base_series_panel(grid, response, cols, labels = paste0("B", seq_len(nb)),
                         main = title, onsets = marks, onset_col = onset_col, ...)
    }
    df <- data.frame(time = grid, response)
    colnames(df)[-1] <- paste0("basis_", seq_len(nb))
  } else {
    .base_series_panel(grid, response, hrf_palette(1), main = title,
                       onsets = marks, onset_col = onset_col, ...)
    df <- data.frame(time = grid, response = response)
  }
  invisible(df)
}

#' Compare Multiple Regressor Objects
#'
#' Plots one or more regressors on a shared time axis. Regressors built from a
#' basis set are expanded into one curve per basis function, and a
#' `regressor_set()` is expanded into one curve per condition. Event onsets
#' (and durations) are drawn as bars under the curves; zero-duration events get
#' a minimum bar width of 0.4% of the time range so that they stay visible. Uses ggplot2 when
#' available, otherwise base graphics. Colours come from [hrf_palette()].
#'
#' Inside a knitr document the result is returned visibly and printed by
#' knitr, which lets document themes (for example dark-mode figure twins)
#' handle the ggplot. At the console the plot is drawn immediately and the
#' data are returned invisibly.
#'
#' @md
#' @param ... Regressor objects (`Reg`) or a `RegSet` to compare. Can be
#'   passed as individual arguments or as a named list.
#' @param grid Numeric vector of time points for evaluation. If NULL (default),
#'   a 0.25 s grid covering all regressors is used. Use a fine grid for the
#'   curve and `samples` to show scan times.
#' @param labels Character vector of labels, one per regressor or one per
#'   plotted curve. If NULL (default), uses list names, condition levels, or
#'   "Regressor_1", "Regressor_2", etc.
#' @param title Character string for the plot title. If NULL (default),
#'   uses "Regressor comparison".
#' @param subtitle Character string for the plot subtitle. If NULL, no subtitle.
#' @param show_onsets Logical or character. If TRUE, mark events for every
#'   regressor in its own colour. If "first", mark only the first regressor's
#'   events, in grey. If FALSE, hide event marks. If NULL (default), uses
#'   "first" for `layout = "overlay"` and TRUE for `layout = "stack"`, so each
#'   panel shows its own events.
#' @param onset_alpha Alpha transparency for event marks. Default is 0.8.
#' @param precision Numeric sampling precision for HRF evaluation. If NULL
#'   (default), the grid spacing capped at 0.33 s.
#' @param use_ggplot Logical; if TRUE and ggplot2 is available, use ggplot2
#'   for plotting. If FALSE, use base R graphics. Default is TRUE.
#' @param draw Logical; draw the plot (default TRUE). With `use_ggplot = TRUE`,
#'   set FALSE to customize the returned data frame's `"plot"` attribute.
#' @param basis Either `"all"` (default) to plot every basis column of a
#'   multi-basis regressor, or `"first"` for the first column only.
#' @param layout Either `"overlay"` (default) or `"stack"` (one panel per
#'   curve on a shared time axis). Stacked panels mark their own events in
#'   grey.
#' @param scales For `layout = "stack"`: `"free_y"` (default) gives each panel
#'   its own y range; `"fixed"` shares one range so amplitudes can be compared.
#' @param samples Optional numeric vector of sample times (for example scan
#'   acquisition times). The regressors are evaluated there and drawn as
#'   points on the curves.
#' @param palette Colour palette: `"auto"` (default), `"categorical"`, or
#'   `"ordered"`; see [plot_hrfs()].
#' @return A data frame in long format with columns 'time', 'Regressor', and
#'   'response'. With ggplot2, the plot is stored in the `"plot"` attribute.
#'   The value is invisible except when drawing inside knitr, where it prints
#'   like a ggplot object (see [plot_hrfs()]).
#' @examples
#' # Create regressors with different HRFs
#' onsets <- c(10, 30, 50)
#' reg1 <- regressor(onsets, HRF_SPMG1)
#' reg2 <- regressor(onsets, HRF_GAMMA)
#' reg3 <- regressor(onsets, HRF_GAUSSIAN)
#'
#' # Compare regressors
#' plot_regressors(reg1, reg2, reg3,
#'                 labels = c("SPM Canonical", "Gamma", "Gaussian"))
#'
#' # Show the scan-time samples of a regressor (TR = 2 s)
#' plot_regressors(reg1, samples = seq(0, 80, by = 2), labels = "SPMG1")
#'
#' # One panel per basis function of a basis-set regressor
#' plot_regressors(regressor(c(10, 40), HRF_SPMG3), layout = "stack")
#'
#' # Compare original vs shifted regressor
#' reg_shifted <- shift(reg1, 5)
#' plot_regressors(reg1, reg_shifted, labels = c("Original", "Shifted +5s"),
#'                 show_onsets = TRUE)
#' @export
plot_regressors <- function(..., grid = NULL, labels = NULL,
                            title = NULL, subtitle = NULL,
                            show_onsets = NULL, onset_alpha = 0.8,
                            precision = NULL, use_ggplot = TRUE, draw = TRUE,
                            basis = c("all", "first"),
                            layout = c("overlay", "stack"), samples = NULL,
                            palette = c("auto", "categorical", "ordered"),
                            scales = c("free_y", "fixed")) {
  scales <- match.arg(scales)
  basis <- match.arg(basis)
  layout <- match.arg(layout)
  palette <- match.arg(palette)
  regs <- list(...)

  if (length(regs) == 1 && inherits(regs[[1]], "RegSet")) {
    set <- regs[[1]]
    regs <- stats::setNames(set$regs, set$levels)
  } else if (length(regs) == 1 && is.list(regs[[1]]) && !inherits(regs[[1]], "Reg")) {
    regs <- regs[[1]]
  }

  n_regs <- length(regs)
  if (n_regs == 0) {
    stop("At least one Reg object must be provided")
  }
  for (i in seq_along(regs)) {
    if (!inherits(regs[[i]], "Reg")) {
      stop("All arguments must be Reg objects. Argument ", i, " is not a Reg.")
    }
  }

  if (is.null(grid)) {
    max_time <- max(vapply(regs, function(r) {
      if (length(r$onsets) == 0) return(as.numeric(r$span))
      max(r$onsets, na.rm = TRUE) + r$span
    }, numeric(1)))
    grid <- seq(0, max_time, by = 0.25)
  }

  if (is.null(show_onsets)) {
    show_onsets <- if (layout == "stack") TRUE else "first"
  }
  precision <- .plot_precision(precision, grid)
  eval_reg <- function(r, t) {
    y <- as.matrix(evaluate(r, t, precision = precision))
    if (basis == "first") y[, 1, drop = FALSE] else y
  }
  mats <- lapply(regs, eval_reg, t = grid)
  nb <- vapply(mats, ncol, integer(1))
  n_series <- sum(nb)
  single_set <- n_regs == 1 && nb[1] > 1

  if (is.null(labels)) {
    base <- names(regs)
    if (is.null(base) || any(!nzchar(base))) base <- paste0("Regressor_", seq_len(n_regs))
  } else if (length(labels) == n_series) {
    base <- NULL
  } else if (length(labels) == n_regs) {
    base <- labels
  } else {
    stop("Length of 'labels' must match number of regressors")
  }
  series_labels <- if (is.null(base)) {
    as.character(labels)
  } else {
    unlist(lapply(seq_len(n_regs), function(i) {
      if (nb[i] == 1) return(base[i])
      if (single_set) paste0("B", seq_len(nb[i])) else paste0(base[i], " B", seq_len(nb[i]))
    }))
  }
  if (any(duplicated(series_labels))) {
    series_labels <- make.unique(series_labels, sep = " ")
  }
  owner <- rep(seq_len(n_regs), nb)

  resp <- do.call(cbind, mats)
  df <- data.frame(
    time = rep(grid, n_series),
    Regressor = factor(rep(series_labels, each = length(grid)), levels = series_labels),
    response = as.vector(resp)
  )

  # Event marks, one row per event and curve.
  onset_data <- NULL
  if (!isFALSE(show_onsets)) {
    first_only <- identical(show_onsets, "first")
    rows <- lapply(seq_len(n_series), function(s) {
      r <- regs[[if (first_only) 1 else owner[s]]]
      if (length(r$onsets) == 0) return(NULL)
      dur <- r$duration %||% 0
      data.frame(onset = r$onsets, duration = rep_len(dur, length(r$onsets)),
                 series = series_labels[s])
    })
    onset_data <- do.call(rbind, rows)
    if (!is.null(onset_data)) {
      # Only mark events that fall inside the plotted window.
      onset_data <- onset_data[onset_data$onset >= min(grid) &
                                 onset_data$onset <= max(grid), , drop = FALSE]
      if (first_only) {
        # Only the first regressor's curve(s) get its event marks.
        onset_data <- onset_data[onset_data$series %in% series_labels[owner == 1], , drop = FALSE]
      }
      attr(onset_data, "by_series") <- !first_only
      attr(onset_data, "alpha") <- onset_alpha
    }
  }

  sample_data <- NULL
  if (!is.null(samples)) {
    smat <- do.call(cbind, lapply(regs, eval_reg, t = samples))
    sample_data <- data.frame(
      time = rep(samples, n_series),
      series = factor(rep(series_labels, each = length(samples)), levels = series_labels),
      response = as.vector(smat)
    )
  }

  if (is.null(title)) {
    title <- if (n_series == 1) series_labels[1] else "Regressor comparison"
  }
  pal <- .series_palette(n_series, palette, ordered_hint = single_set && nb[1] > 3)

  if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    gdf <- df
    names(gdf)[2] <- "series"
    p <- .gg_series(gdf, pal$colours, title = title, subtitle = subtitle,
                    layout = layout, onsets = onset_data, samples = sample_data,
                    linetypes = pal$type == "categorical" && n_series > 4,
                    scales = scales)
    return(.finish_gg(df, p, draw))
  }
  if (draw) {
    marks <- if (!is.null(onset_data)) unique(onset_data$onset) else NULL
    .base_series_panel(grid, resp, pal$colours, labels = series_labels,
                       main = title, onsets = marks,
                       onset_col = grDevices::adjustcolor(.hrf_neutral, alpha.f = onset_alpha))
    if (!is.null(sample_data)) {
      graphics::points(sample_data$time, sample_data$response, pch = 19, cex = 0.6,
                       col = pal$colours[as.integer(sample_data$series)])
    }
  }
  invisible(df)
}

#' Plot a Feature Regressor
#'
#' @md
#' @param x A `FeatureReg` object created by [feature_regressor()].
#' @param grid Numeric vector of time points for evaluation. If `NULL`
#'   (default), a grid from 0 to max(times) + span with step 0.25 s is used.
#' @param show_onsets Logical; if `TRUE`, mark sample times with ticks on the
#'   time axis. Defaults to `FALSE` because a dense feature has one sample per
#'   bin.
#' @param onset_color Colour for sample-time ticks. If NULL (default), a
#'   neutral grey.
#' @param onset_alpha Alpha transparency for sample-time ticks. Default is 0.5.
#' @param precision Numeric sampling precision for HRF evaluation. If NULL
#'   (default), the grid spacing capped at 0.33 s.
#' @param layout For multi-basis HRFs, `"stack"` (default) or `"overlay"`.
#' @param ... Additional arguments passed to the underlying plot functions.
#' @return Invisibly returns a data frame with the time and response values.
#' @examples
#' feat <- feature_regressor(abs(sin(seq(0, 8, by = 0.1))), dt = 0.1)
#' plot(feat, grid = seq(0, 12, by = 0.5))
#' @method plot FeatureReg
#' @export
plot.FeatureReg <- function(x, grid = NULL, show_onsets = FALSE,
                            onset_color = NULL, onset_alpha = 0.5,
                            precision = NULL, layout = c("stack", "overlay"), ...) {
  layout <- match.arg(layout)
  if (is.null(grid)) {
    max_time <- max(x$onsets, na.rm = TRUE) + x$span
    grid <- seq(0, max_time, by = 0.25)
  }
  response <- evaluate(x, grid, precision = .plot_precision(precision, grid))
  hrf_name <- attr(x$hrf, "name") %||% "custom HRF"
  title <- sprintf("Feature regressor (%s), %d samples", hrf_name, length(x$onsets))
  .plot_reg_base(response, grid, x$onsets, title, show_onsets, onset_color,
                 onset_alpha, layout, ...)
}

# Evaluation precision for plotting: the grid spacing, capped at 0.33 s (the
# evaluate() default) and floored at 0.01 s.
.plot_precision <- function(precision, grid) {
  if (!is.null(precision)) return(precision)
  # signif() removes floating-point residue from seq() spacing (0.0499999...)
  step <- if (length(grid) > 1) signif(min(diff(sort(unique(grid)))), 8) else 0.33
  min(0.33, max(0.01, step))
}

# Colours for the columns of one basis set: categorical for up to three
# functions (e.g. canonical + derivatives), ordered beyond that.
.basis_palette <- function(nb) {
  hrf_palette(nb, if (nb > 3) "ordered" else "categorical")
}
