## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, 
  fig.height = 5,
  message = FALSE,
  warning = FALSE
)
library(fmrihrf)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----basic_regressor----------------------------------------------------------
# Define event onsets
onsets <- seq(0, 10 * 12, by = 12)

# Create the regressor object
# Uses HRF_SPMG1 by default if no hrf is specified
# Duration is 0 by default
reg1 <- regressor(onsets = onsets, hrf = HRF_SPMG1)

# Print the regressor object to see its properties (uses new print.Reg method)
print(reg1)

# Access components using helper functions
head(onsets(reg1))
nbasis(reg1)

## ----evaluate_plot_basic------------------------------------------------------
# Define a time grid corresponding to scan times (e.g., TR=2s)
TR <- 2
scan_times <- seq(0, 140, by = TR)

# Plot the regressor using the plot() method
# This automatically evaluates and plots with onset markers
plot(reg1, grid = scan_times)

## ----varying_duration---------------------------------------------------------
# Example onsets and durations
onsets_var_dur <- seq(0, 5 * 12, length.out = 6)
durations_var <- 1:length(onsets_var_dur) # Durations increase from 1s to 6s

# Create regressor with varying durations
reg_var_dur <- regressor(onsets_var_dur, HRF_SPMG1, duration = durations_var)
print(reg_var_dur)

# Plot the regressor
scan_times_dur <- seq(0, max(onsets_var_dur) + 30, by = TR)
plot(reg_var_dur, grid = scan_times_dur)

## ----duration_no_summate------------------------------------------------------
# Create regressor with varying durations, summate=FALSE
reg_var_dur_nosum <- regressor(onsets_var_dur, HRF_SPMG1,
                               duration = durations_var, summate = FALSE)

# Compare summating vs non-summating using plot_regressors()
plot_regressors(reg_var_dur, reg_var_dur_nosum,
                labels = c("summate=TRUE", "summate=FALSE"),
                grid = scan_times_dur,
                title = "Effect of Summation on Response",
                subtitle = "Same events with varying durations")

## ----parametric_modulation----------------------------------------------------
# Example onsets and amplitudes (e.g., representing task difficulty)
onsets_amp <- seq(0, 10 * 12, length.out = 11)
amplitudes_raw <- 1:length(onsets_amp)

# It's common practice to center the modulator
amplitudes_scaled <- scale(amplitudes_raw, center = TRUE, scale = FALSE)

# Create the parametric regressor
reg_amp <- regressor(onsets_amp, HRF_SPMG1, amplitude = amplitudes_scaled)
print(reg_amp)

# Plot the parametric regressor
scan_times_amp <- seq(0, max(onsets_amp) + 30, by = TR)
plot(reg_amp, grid = scan_times_amp)

## ----duration_amplitude-------------------------------------------------------
set.seed(123)
onsets_comb <- seq(0, 10 * 12, length.out = 11)
amps_comb <- scale(1:length(onsets_comb), center = TRUE, scale = FALSE)
durs_comb <- sample(1:5, length(onsets_comb), replace = TRUE)

reg_comb <- regressor(onsets_comb, HRF_SPMG1, 
                      amplitude = amps_comb, duration = durs_comb)
print(reg_comb)

# Plot the combined regressor
scan_times_comb <- seq(0, max(onsets_comb) + 30, by = TR)
plot(reg_comb, grid = scan_times_comb)

## ----basis_set_regressor------------------------------------------------------
# Use a B-spline basis set
onsets_basis <- seq(0, 10 * 12, length.out = 11)
hrf_basis <- HRF_BSPLINE # Uses N=5 basis functions by default

reg_basis <- regressor(onsets_basis, hrf_basis)
print(reg_basis)
nbasis(reg_basis) # Should be 5

# Evaluate - this returns a matrix
scan_times_basis <- seq(0, max(onsets_basis) + 30, by = TR)
pred_basis_matrix <- evaluate(reg_basis, scan_times_basis)
dim(pred_basis_matrix) # rows = time points, cols = basis functions

# Plot using the plot() method - automatically handles multi-basis
plot(reg_basis, grid = scan_times_basis)

## ----shift_regressor----------------------------------------------------------
# Original regressor
reg_orig <- regressor(onsets = c(10, 30, 50), hrf = HRF_SPMG1)

# Shifted regressor (delay by 5 seconds)
reg_shifted <- shift(reg_orig, shift_amount = 5)

onsets(reg_orig)
onsets(reg_shifted) # Onsets are now 15, 35, 55

# Compare original and shifted using plot_regressors()
scan_times_shift <- seq(0, 80, by = TR)
plot_regressors(reg_orig, reg_shifted,
                labels = c("Original", "Shifted +5s"),
                grid = scan_times_shift,
                show_onsets = TRUE,  # Show onsets for both
                title = "Shifting a Regressor")

