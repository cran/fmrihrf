## ----setup, include = FALSE--------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  # Sharper website figures; keep the CRAN vignette archive compact.
  fig.retina = if (identical(Sys.getenv("IN_PKGDOWN"), "true")) 2 else 1,
  fig.width = 7,
  fig.height = 4,
  message = FALSE,
  warning = FALSE
)
# CRAN builds: skip dark-mode figure twins to keep the source package small;
# the pkgdown site (IN_PKGDOWN = "true") keeps them.
if (!identical(Sys.getenv("IN_PKGDOWN"), "true")) {
  options(albersdown.dark_figures = FALSE)
}
library(fmrihrf)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----basic_regressor---------------------------------------------
# Define event onsets
onsets <- seq(0, 10 * 12, by = 12)

# Create the regressor object
# Uses HRF_SPMG1 by default if no hrf is specified
# Duration is 0 by default
reg1 <- regressor(onsets = onsets, hrf = HRF_SPMG1)

# Access components using helper functions
head(onsets(reg1))
nbasis(reg1)

## ----convolution_idea, echo = FALSE, fig.height = 4.6, fig.alt="Three stacked panels. Top: eleven unit event sticks every 12 seconds. Middle: the SPM canonical HRF over 24 seconds. Bottom: the regressor, the sum of one HRF copy per event, which rises and falls after each event and overlaps slightly with the next."----
conv_grid <- seq(0, 150, by = 0.1)
hrf_grid <- seq(0, 24, by = 0.1)
conv_df <- rbind(
  data.frame(time = rep(onsets, each = 3), value = rep(c(0, 1, NA), length(onsets)),
             panel = "1. Events (onsets)"),
  data.frame(time = hrf_grid, value = HRF_SPMG1(hrf_grid) / max(HRF_SPMG1(hrf_grid)),
             panel = "2. HRF for one event (peak = 1)"),
  data.frame(time = conv_grid, value = evaluate(reg1, conv_grid) / max(HRF_SPMG1(hrf_grid)),
             panel = "3. Regressor = events convolved with HRF")
)
conv_df$panel <- factor(conv_df$panel, levels = unique(conv_df$panel))
ggplot(conv_df, aes(time, value, colour = panel)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(linetype = panel), linewidth = 0.9, na.rm = TRUE) +
  scale_linetype_manual(values = c("solid", "22", "solid")) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("grey45", "grey45", hrf_palette(1))) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(title = "From events to a regressor", x = "Time (s)", y = NULL) +
  theme(legend.position = "none", strip.text = element_text(hjust = 0),
        plot.title.position = "plot")

## ----evaluate_basic----------------------------------------------
# Define a time grid corresponding to scan times (e.g., TR=2s)
TR <- 2
scan_times <- seq(0, 140, by = TR)

# One value per scan
head(evaluate(reg1, scan_times))

## ----evaluate_plot_basic, fig.alt="Modeled SPM response to eleven events spaced 12 seconds apart, drawn as a smooth curve with dots at the 2 second scan times. Grey bars below the curve mark the event onsets.", fig.cap="The curve is the modeled response; dots are its values at the scan times (TR = 2 s); grey bars mark event onsets."----
plot_regressors(reg1, grid = seq(0, 140, by = 0.1), samples = scan_times,
                labels = "SPMG1 regressor",
                title = "Predicted response to 11 events",
                subtitle = "Evaluated at every scan (TR = 2 s)")

## ----varying_duration--------------------------------------------
# Example onsets and durations
onsets_var_dur <- seq(0, 5 * 12, length.out = 6)
durations_var <- 1:length(onsets_var_dur) # Durations increase from 1s to 6s

# Create regressor with varying durations
reg_var_dur <- regressor(onsets_var_dur, HRF_SPMG1, duration = durations_var)

scan_times_dur <- seq(0, max(onsets_var_dur) + 30, by = TR)
fine_grid_dur <- seq(0, max(scan_times_dur), by = 0.1)

## ----varying_duration_plot, fig.alt="Regressor for six events whose durations grow from 1 to 6 seconds. Bars below the curve show each event's duration; response peaks grow with duration."----
plot_regressors(reg_var_dur, grid = fine_grid_dur,
                labels = "Durations 1-6 s",
                title = "Increasing event durations")

## ----duration_no_summate, fig.alt="Summed and averaged regressors for six events with durations 1 to 6 seconds. The summed response grows with each longer event; the averaged response stays at about the same height."----
# Create regressor with varying durations, summate=FALSE
reg_var_dur_nosum <- regressor(onsets_var_dur, HRF_SPMG1,
                               duration = durations_var, summate = FALSE)

# Compare summating vs non-summating using plot_regressors()
plot_regressors(reg_var_dur, reg_var_dur_nosum,
                labels = c("summate = TRUE", "summate = FALSE"),
                grid = fine_grid_dur,
                title = "Summed vs. averaged responses",
                subtitle = "Same six events, durations 1-6 s")

## ----parametric_modulation---------------------------------------
# Example onsets and amplitudes (e.g., representing task difficulty)
onsets_amp <- seq(0, 10 * 12, length.out = 11)
amplitudes_raw <- 1:length(onsets_amp)

# It's common practice to center the modulator
amplitudes_scaled <- scale(amplitudes_raw, center = TRUE, scale = FALSE)

# Create the parametric regressor
reg_amp <- regressor(onsets_amp, HRF_SPMG1, amplitude = amplitudes_scaled)

# The centred amplitudes run from -5 to 5; the middle event has amplitude 0
drop(amplitudes_scaled)

## ----parametric_modulation_plot, fig.alt="Parametric regressor for 11 events with centred amplitudes from -5 to 5. The response dips below zero for early events and rises above zero for late events; the middle event produces no response."----
fine_grid_amp <- seq(0, max(onsets_amp) + 30, by = 0.1)
plot_regressors(reg_amp, grid = fine_grid_amp,
                labels = "Amplitude-modulated",
                title = "Parametric modulation",
                subtitle = "Mean-centred amplitudes from -5 to 5")

## ----feature_regressor-------------------------------------------
dt <- 0.1
feat_times <- seq(0, 20, by = dt)
# Simulated acoustic envelope
rms <- abs(sin(2 * pi * feat_times / 8)) * (0.5 + 0.5 * sin(2 * pi * feat_times / 20))

feat <- feature_regressor(rms, dt = dt, hrf = HRF_SPMG1)

## ----feature_regressor_plot, fig.height = 4.6, fig.alt="Two panels. Top: the centred acoustic envelope, oscillating around zero over 20 seconds. Bottom: the predicted BOLD response, a smoothed and delayed version that continues for about 20 seconds after the feature ends."----
feat_grid <- seq(0, max(feat_times) + 30, by = 0.1)
feat_df <- rbind(
  data.frame(time = feat_times, value = rms - mean(rms),
             panel = "Feature (centred), input"),
  data.frame(time = feat_grid, value = evaluate(feat, feat_grid, precision = dt),
             panel = "Predicted BOLD, output")
)
feat_df$panel <- factor(feat_df$panel, levels = unique(feat_df$panel))
ggplot(feat_df, aes(time, value, colour = panel)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("grey45", hrf_palette(1))) +
  scale_y_continuous(breaks = function(l) {
    b <- pretty(l, n = 3)
    b[b >= l[1] & b <= l[2]]
  }) +
  labs(title = "Feature regressor", subtitle = "Input (centred feature) and output",
       x = "Time (s)", y = NULL) +
  theme(legend.position = "none", strip.text = element_text(hjust = 0),
        plot.title.position = "plot")

## ----duration_amplitude, fig.alt="Regressor for 11 events with random durations of 1 to 5 seconds and centred amplitudes from -5 to 5. Bars below the curve show the durations; responses go from negative to positive over the run."----
set.seed(123)
onsets_comb <- seq(0, 10 * 12, length.out = 11)
amps_comb <- scale(1:length(onsets_comb), center = TRUE, scale = FALSE)
durs_comb <- sample(1:5, length(onsets_comb), replace = TRUE)

reg_comb <- regressor(onsets_comb, HRF_SPMG1, 
                      amplitude = amps_comb, duration = durs_comb)

fine_grid_comb <- seq(0, max(onsets_comb) + 30, by = 0.1)
plot_regressors(reg_comb, grid = fine_grid_comb,
                labels = "Duration and amplitude",
                title = "Duration and amplitude modulation",
                subtitle = "Bar width = duration")

## ----basis_set_regressor-----------------------------------------
# Use a B-spline basis set
onsets_basis <- seq(0, 10 * 12, length.out = 11)
hrf_basis <- HRF_BSPLINE # Uses N=5 basis functions by default

reg_basis <- regressor(onsets_basis, hrf_basis)
nbasis(reg_basis) # Should be 5

# Evaluate - this returns a matrix
scan_times_basis <- seq(0, max(onsets_basis) + 30, by = TR)
pred_basis_matrix <- evaluate(reg_basis, scan_times_basis)
dim(pred_basis_matrix) # rows = time points, cols = basis functions

## ----basis_set_regressor_plot, fig.height = 5.5, fig.alt="Five stacked panels, one per B-spline basis function, showing the regressor columns over the first 58 seconds. B1 peaks shortly after each event and each later basis function peaks later."----
plot_regressors(reg_basis, grid = seq(0, 58, by = 0.1), layout = "stack",
                title = "One column per basis function",
                subtitle = "B-spline basis (N = 5), events every 12 s")

## ----shift_regressor---------------------------------------------
# Original regressor
reg_orig <- regressor(onsets = c(10, 30, 50), hrf = HRF_SPMG1)

# Shifted regressor (delay by 5 seconds)
reg_shifted <- shift(reg_orig, shift_amount = 5)

onsets(reg_orig)
onsets(reg_shifted) # Onsets are now 15, 35, 55


## ----shift_regressor_plot, fig.alt="Original and shifted regressors for three events. The shifted curve and its onset bars are 5 seconds later; the peak heights are identical."----
plot_regressors(reg_orig, reg_shifted,
                labels = c("Original", "Shifted +5 s"),
                grid = seq(0, 80, by = 0.1),
                show_onsets = TRUE,  # Show onsets for both
                title = "Shifting a regressor")

