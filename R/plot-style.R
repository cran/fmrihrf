#' Colour palettes for HRF and regressor plots
#'
#' `hrf_palette()` returns the colours used by every plotting function in
#' fmrihrf, so figures made with [plot_hrfs()], [plot_regressors()], and the
#' `plot()` methods share one visual language.
#'
#' Two palette types are provided:
#'
#' * `"categorical"`: up to six distinct hues for unordered series (conditions,
#'   HRF families). Every colour has at least 3:1 contrast against both white
#'   and near-black backgrounds, and the first four stay distinguishable under
#'   simulated deuteranopia, protanopia, and tritanopia. When more than six
#'   colours are requested the ordered palette is used instead.
#' * `"ordered"`: a hue ramp (violet, blue, teal, green, ochre) for series
#'   with a natural order, such as the functions of a basis set, lags, or
#'   event durations. Lightness is held in the same mid range, so the ramp is
#'   readable on light and dark backgrounds. Brick red, the first
#'   categorical colour, is left out of the ramp so that it keeps one meaning.
#'
#' @md
#' @param n Number of colours. If `NULL`, returns the six categorical colours
#'   (or five ramp anchors for `type = "ordered"`).
#' @param type Either `"categorical"` or `"ordered"`.
#' @return A character vector of hex colours.
#' @seealso [scale_colour_hrf()] for the matching ggplot2 scales.
#' @examples
#' hrf_palette()
#' hrf_palette(3)
#' hrf_palette(8, type = "ordered")
#' @export
hrf_palette <- function(n = NULL, type = c("categorical", "ordered")) {
  type <- match.arg(type)
  categorical <- c(
    brick = "#BD4B3F", cerulean = "#1A95AE", violet = "#7A51C8",
    ochre = "#B08214", green = "#2B8667", mauve = "#9C6687"
  )
  anchors <- c("#7A51C8", "#3F6FD0", "#1A95AE", "#2B8667", "#B08214")
  if (is.null(n)) {
    return(if (type == "categorical") unname(categorical) else anchors)
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("'n' must be a single non-negative integer")
  }
  if (n == 0L) {
    return(character(0))
  }
  if (type == "categorical" && n <= length(categorical)) {
    return(unname(categorical[seq_len(n)]))
  }
  if (n == 1L) {
    return(unname(categorical[1]))
  }
  grDevices::colorRampPalette(anchors, space = "Lab")(n)
}

#' ggplot2 colour scales matching fmrihrf plots
#'
#' Discrete colour and fill scales built on [hrf_palette()]. Use them in your
#' own ggplot2 figures to match the output of [plot_hrfs()] and
#' [plot_regressors()].
#'
#' @md
#' @param type Either `"categorical"` or `"ordered"`; see [hrf_palette()].
#' @param ... Further arguments passed to [ggplot2::discrete_scale()], such as
#'   `name`, `labels`, or `guide`.
#' @param aesthetics The aesthetics the scale applies to.
#' @return A ggplot2 scale object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   t <- seq(0, 24, by = 0.2)
#'   df <- data.frame(
#'     time = rep(t, 3),
#'     response = c(HRF_SPMG1(t), HRF_GAMMA(t), HRF_GAUSSIAN(t)),
#'     hrf = rep(c("SPMG1", "Gamma", "Gaussian"), each = length(t))
#'   )
#'   ggplot2::ggplot(df, ggplot2::aes(time, response, colour = hrf)) +
#'     ggplot2::geom_line() +
#'     scale_colour_hrf()
#' }
#' @export
scale_colour_hrf <- function(type = c("categorical", "ordered"), ...,
                             aesthetics = "colour") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("scale_colour_hrf() requires the 'ggplot2' package")
  }
  type <- match.arg(type)
  ggplot2::discrete_scale(aesthetics, palette = function(n) hrf_palette(n, type), ...)
}

#' @md
#' @rdname scale_colour_hrf
#' @export
scale_color_hrf <- scale_colour_hrf

#' @md
#' @rdname scale_colour_hrf
#' @export
scale_fill_hrf <- function(type = c("categorical", "ordered"), ...,
                           aesthetics = "fill") {
  scale_colour_hrf(type = type, ..., aesthetics = aesthetics)
}

# ---------------------------------------------------------------------------
# Internal plotting helpers shared by plot_hrfs(), plot_regressors() and the
# base-graphics plot() methods.
# ---------------------------------------------------------------------------

# Neutral ink for zero lines, onset marks and reference curves. albersdown's
# dark figure twins mirror these greys automatically.
.hrf_neutral <- "grey45"
.hrf_rule <- "grey70"

# Wrap text that would overflow a phone-width figure into balanced lines, so
# desktop renders do not end with a single orphaned word.
.wrap_title <- function(x, width = 38) {
  if (is.null(x)) {
    return(NULL)
  }
  if (nchar(x) <= width) {
    return(x)
  }
  n_lines <- ceiling(nchar(x) / width)
  target <- min(width, ceiling(nchar(x) / n_lines) + 4)
  paste(strwrap(x, target), collapse = "\n")
}

# Resolve the palette type for a set of series.
.series_palette <- function(n, palette = c("auto", "categorical", "ordered"),
                            ordered_hint = FALSE) {
  palette <- match.arg(palette)
  if (palette == "auto") {
    palette <- if (ordered_hint || n > 6) "ordered" else "categorical"
  }
  list(type = palette, colours = hrf_palette(n, palette))
}

# Wrap the long-form data returned by the comparison helpers. Inside knitr the
# ggplot is not printed here; the value is returned visibly so knitr's
# knit_print() machinery (and themes that add dark-mode twins) handles it.
.finish_gg <- function(df, p, draw) {
  attr(df, "plot") <- p
  if (!draw) {
    return(invisible(df))
  }
  if (isTRUE(getOption("knitr.in.progress"))) {
    class(df) <- c("fmrihrf_plot", class(df))
    return(df)
  }
  print(p)
  invisible(df)
}

#' @md
#' @export
`[.fmrihrf_plot` <- function(x, ...) {
  .plain_df(x)[...]
}

#' @md
#' @export
`$<-.fmrihrf_plot` <- function(x, name, value) {
  x <- .plain_df(x)
  x[[name]] <- value
  x
}

#' @md
#' @export
`[[<-.fmrihrf_plot` <- function(x, ..., value) {
  x <- .plain_df(x)
  x[[...]] <- value
  x
}

#' @md
#' @export
`names<-.fmrihrf_plot` <- function(x, value) {
  x <- .plain_df(x)
  names(x) <- value
  x
}

#' @export
`[<-.fmrihrf_plot` <- function(x, ..., value) {
  x <- .plain_df(x)
  x[...] <- value
  x
}

#' @md
#' @export
as.data.frame.fmrihrf_plot <- function(x, ...) {
  .plain_df(x)
}

#' @md
#' @export
print.fmrihrf_plot <- function(x, ...) {
  p <- attr(x, "plot")
  if (inherits(p, "ggplot")) {
    print(p)
    return(invisible(x))
  }
  class(x) <- setdiff(class(x), "fmrihrf_plot")
  print(x, ...)
  invisible(x)
}

# Registered for knitr::knit_print in .onLoad() when knitr is available.
knit_print.fmrihrf_plot <- function(x, ...) {
  p <- attr(x, "plot")
  if (inherits(p, "ggplot")) {
    return(knitr::knit_print(p, ...))
  }
  class(x) <- setdiff(class(x), "fmrihrf_plot")
  knitr::knit_print(x, ...)
}

# Legend columns that fit a phone-width figure (about 36 characters a row).
.legend_ncol <- function(labels) {
  n <- length(labels)
  width <- nchar(labels) + 4
  if (sum(width) <= 34) return(n)
  # Otherwise as many columns of the widest entry as fit in one row.
  max(1L, min(n, floor(34 / max(width))))
}

# y breaks under a label budget, so tick labels stay apart in short (phone)
# panels. Breaks always come from one evenly spaced grid of a nice step, so
# axes never look non-linear. Among grids from fine to coarse, the first one that
#   - has at most `budget` labels (4 for overlays, or 5 when only that labels
#     both signs; 3 for stacked panels whose data cross zero, else 2),
#   - leaves at least `gap` of the range between labels (0.12 overlay,
#     0.22 stacked),
#   - and, when the data cross zero, labels both signs
# is used; if no grid labels both signs, the first grid meeting the other two
# conditions is used. `floor` reserves the bottom fraction of the range for
# event marks, where no labels are placed.
.budget_breaks <- function(limits, stack = FALSE, floor = 0) {
  lower <- limits[1] + floor * diff(limits)
  upper <- limits[2]
  range <- upper - lower
  if (!is.finite(range) || range <= 0) return(pretty(limits))
  # Data must reach meaningfully below zero (10% of the range, beyond the
  # scale's own expansion) to count as crossing it.
  crosses <- lower < -0.1 * range && upper > 0
  budget <- if (stack) (if (crosses) 3 else 2) else 4
  gap <- if (stack) 0.22 else 0.12
  # Candidate grids: multiples of a nice step (1, 2, 2.5 or 5 x 10^k), which
  # include zero whenever the data cross it.
  k <- floor(log10(range))
  steps <- sort(as.vector(outer(c(1, 2, 2.5, 5), 10^((k - 2):(k + 1)))))
  grid <- function(step) {
    b <- seq(ceiling(lower / step - 1e-9), floor(upper / step + 1e-9)) * step
    b[abs(b) < step * 1e-9] <- 0
    b
  }
  fallback <- NULL
  for (step in steps) {
    b <- grid(step)
    if (!length(b)) next
    if (length(b) > budget || (length(b) > 1 && step < gap * range)) next
    if (!crosses || (any(b < 0) && any(b > 0))) return(b)
    if (is.null(fallback)) fallback <- b
  }
  # Overlays may use one extra label when that is the only way to label both
  # sides of zero on an even grid.
  if (crosses && !stack) {
    for (step in steps) {
      b <- grid(step)
      if (length(b) <= budget + 1L && step >= gap * range && any(b < 0) && any(b > 0)) {
        return(b)
      }
    }
  }
  fallback %||% pretty(c(lower, upper), n = 1)
}

# Removing the knitr print class: any modification returns a plain data frame.
.plain_df <- function(x) {
  class(x) <- setdiff(class(x), "fmrihrf_plot")
  attr(x, "plot") <- NULL
  x
}

# Build the shared ggplot for long-form series data.
#   df:        data.frame(time, series (factor), response)
#   onsets:    NULL or data.frame(onset, duration, series (factor or NA))
#   samples:   NULL or data.frame(time, series, response)
.gg_series <- function(df, colours, title = NULL, subtitle = NULL,
                       xlab = "Time (s)", ylab = "Response",
                       layout = c("overlay", "stack"), onsets = NULL,
                       samples = NULL, direct_labels = FALSE,
                       linetypes = FALSE, reference = NULL,
                       scales = "free_y") {
  layout <- match.arg(layout)
  gg <- asNamespace("ggplot2")
  n <- nlevels(df$series)
  names(colours) <- levels(df$series)

  yr <- range(df$response, reference$response, 0, finite = TRUE)
  span_y <- diff(yr)
  if (!is.finite(span_y) || span_y == 0) span_y <- 1
  xr <- range(df$time, finite = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$response,
                                        colour = .data$series)) +
    ggplot2::geom_hline(yintercept = 0, colour = .hrf_rule, linewidth = 0.3)

  # Event marks: a bar from onset to onset + duration, drawn just below the
  # lowest point of the curve(s) it belongs to (per panel when stacked).
  if (!is.null(onsets) && nrow(onsets) > 0) {
    # One set of marks per panel when stacked, so colour is not needed there.
    by_series <- !isFALSE(attr(onsets, "by_series")) && layout == "overlay"
    onsets$series <- factor(onsets$series, levels = levels(df$series))
    width <- pmax(onsets$duration, diff(xr) * 0.004)
    onsets$xend <- onsets$onset + width
    if (layout == "stack" && scales == "fixed") {
      # Shared y axis: place every panel's marks from the shared range so the
      # label floor below holds for all panels.
      onsets$y0 <- yr[1] - 0.22 * span_y
      onsets$y1 <- yr[1] - 0.06 * span_y
    } else if (layout == "stack") {
      lo <- tapply(df$response, df$series, function(v) min(v, 0, na.rm = TRUE))
      rg <- tapply(df$response, df$series, function(v) {
        r <- diff(range(v, 0, na.rm = TRUE))
        if (!is.finite(r) || r == 0) 1 else r
      })
      key <- as.character(onsets$series)
      onsets$y0 <- lo[key] - 0.22 * rg[key]
      onsets$y1 <- lo[key] - 0.06 * rg[key]
    } else {
      onsets$y0 <- yr[1] - 0.16 * span_y
      onsets$y1 <- yr[1] - 0.06 * span_y
    }
    if (by_series) {
      p <- p + ggplot2::geom_rect(
        data = onsets, inherit.aes = FALSE,
        ggplot2::aes(xmin = .data$onset, xmax = .data$xend,
                     ymin = .data$y0, ymax = .data$y1, fill = .data$series),
        colour = NA, alpha = attr(onsets, "alpha") %||% 1, show.legend = FALSE
      ) + ggplot2::scale_fill_manual(values = colours, guide = "none")
    } else {
      p <- p + ggplot2::geom_rect(
        data = onsets, inherit.aes = FALSE,
        ggplot2::aes(xmin = .data$onset, xmax = .data$xend,
                     ymin = .data$y0, ymax = .data$y1),
        fill = "grey35", colour = NA, alpha = attr(onsets, "alpha") %||% 1
      )
    }
  }

  if (!is.null(reference) && nrow(reference) > 0) {
    p <- p + ggplot2::geom_line(data = reference, inherit.aes = FALSE,
                                ggplot2::aes(x = .data$time, y = .data$response),
                                colour = .hrf_neutral, linewidth = 0.7, linetype = "22")
  }

  if (linetypes && layout == "overlay") {
    p <- p + ggplot2::geom_line(ggplot2::aes(linetype = .data$series), linewidth = 0.9)
  } else {
    p <- p + ggplot2::geom_line(linewidth = 0.9)
  }

  if (!is.null(samples) && nrow(samples) > 0) {
    samples$series <- factor(samples$series, levels = levels(df$series))
    p <- p + ggplot2::geom_point(data = samples, size = 2, stroke = 0)
  }

  p <- p + ggplot2::scale_colour_manual(values = colours, name = NULL)
  if (linetypes && layout == "overlay") {
    p <- p + ggplot2::scale_linetype_discrete(name = NULL)
  }

  if (direct_labels && layout == "overlay") {
    peaks <- do.call(rbind, lapply(split(df, df$series), function(d) {
      if (!nrow(d) || all(!is.finite(d$response))) return(NULL)
      # Centre of the peak (plateaus such as FIR bins get a centred label)
      top <- which(d$response >= max(d$response, na.rm = TRUE) - 1e-6 * diff(yr))
      d[top[ceiling(length(top) / 2)], , drop = FALSE]
    }))
    p <- p + ggplot2::geom_text(
      data = peaks, ggplot2::aes(label = .data$series),
      vjust = -0.6, size = 3.2, fontface = "bold", show.legend = FALSE
    )
  }
  # Fraction of the y range taken by event marks (no tick labels there).
  # Tick labels must stay above the event marks. Marks run from 0.22 (stacked)
  # or 0.16 (overlay) of the data range r below the data down to 0.06 r below
  # it; the scale then adds 8% (stacked) or 5% (overlay) expansion at each
  # end. Top of the marks as a fraction of the axis passed to the breaks
  # function:
  #   stacked: (0.08 * 1.22 + 0.16) / (1.22 * 1.16) = 0.182
  #   overlay: (0.05 * 1.16 + 0.10) / (1.16 * 1.10) = 0.124
  # (with direct labels the overlay top expansion is larger, lowering it).
  has_marks <- !is.null(onsets) && nrow(onsets) > 0
  mark_floor <- if (has_marks) 0.13 else 0
  stack_floor <- if (has_marks) 0.19 else 0
  if (layout == "overlay") {
    # At most four y labels keep tick labels apart in short phone renders.
    p <- p + ggplot2::scale_y_continuous(
      breaks = function(l) .budget_breaks(l, floor = mark_floor),
      expand = ggplot2::expansion(mult = c(0.05, if (direct_labels) 0.14 else 0.05))
    )
  }

  if (layout == "stack") {
    # Compact small multiples: few y breaks and thin strips, so panels stay
    # legible when a theme renders the figure at phone size.
    p <- p + ggplot2::facet_wrap(~series, ncol = 1, scales = scales,
                                 strip.position = "top") +
      ggplot2::scale_y_continuous(breaks = function(l) .budget_breaks(l, stack = TRUE, floor = stack_floor),
                                  expand = ggplot2::expansion(mult = 0.08)) +
      ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0, size = ggplot2::rel(0.85),
                                                        margin = ggplot2::margin(1, 2, 1, 2)),
                     legend.position = "none",
                     panel.spacing.y = grid::unit(0.35, "lines"))
  } else if (direct_labels || n == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    ncol <- .legend_ncol(levels(df$series))
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = ncol, byrow = TRUE),
      linetype = if (linetypes) ggplot2::guide_legend(ncol = ncol, byrow = TRUE) else "none"
    ) + ggplot2::theme(legend.position = "bottom", legend.justification = "left")
  }

  p + ggplot2::labs(title = .wrap_title(title), subtitle = .wrap_title(subtitle, 42),
                    x = xlab, y = ylab) +
    ggplot2::theme(plot.title.position = "plot")
}

# Base-graphics counterpart: draw one panel of series with the shared palette.
.base_series_panel <- function(time, Y, colours, labels = NULL, main = NULL,
                               xlab = "Time (s)", ylab = "Response",
                               onsets = NULL, onset_col = .hrf_neutral,
                               legend = TRUE, ...) {
  Y <- as.matrix(Y)
  args <- list(...)
  ylim <- args$ylim %||% range(Y, 0, finite = TRUE)
  args$ylim <- NULL
  do.call(graphics::plot, c(list(x = range(time), y = ylim, type = "n",
                                 xlab = xlab, ylab = ylab, main = ""), args))
  graphics::abline(h = 0, col = .hrf_rule, lwd = 0.8)
  if (!is.null(onsets) && length(onsets) > 0) {
    graphics::rug(onsets, ticksize = 0.035, side = 1, col = onset_col, lwd = 1.5)
  }
  for (j in seq_len(ncol(Y))) {
    graphics::lines(time, Y[, j], col = colours[j], lwd = 1.8)
  }
  if (!is.null(main)) {
    graphics::title(main = main, adj = 0, line = if (legend && ncol(Y) > 1) 2.2 else 0.8,
                    cex.main = 1, font.main = 2)
  }
  if (legend && ncol(Y) > 1 && !is.null(labels)) {
    graphics::legend("bottomleft", legend = labels, col = colours, lwd = 1.8,
                     bty = "n", horiz = TRUE, inset = c(0, 1), xpd = TRUE,
                     cex = 0.8, seg.len = 1.2, x.intersp = 0.5)
  }
  invisible(NULL)
}
