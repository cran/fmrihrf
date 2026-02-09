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
library(dplyr) # For pipe operator %>%
library(ggplot2) # For plotting
library(tidyr) # For data manipulation

## ----print_hrfs---------------------------------------------------------------
# SPM canonical HRF (based on difference of two gamma functions)
print(HRF_SPMG1)

# Gaussian HRF
print(HRF_GAUSSIAN)

## ----evaluate_basic_hrfs------------------------------------------------------
time_points <- seq(0, 25, by = 0.1)

# Compare HRFs using plot_hrfs() - normalize = TRUE scales to peak at 1.0
plot_hrfs(HRF_SPMG1, HRF_GAUSSIAN,
          labels = c("SPM Canonical", "Gaussian"),
          normalize = TRUE,
          title = "Comparison of SPM Canonical and Gaussian HRFs",
          subtitle = "HRFs normalized to peak at 1.0 for shape comparison")

## ----modify_gaussian_params---------------------------------------------------
# Create Gaussian HRFs with different parameters using gen_hrf
# Note: hrf_gaussian is the underlying function, not the HRF object HRF_GAUSSIAN
hrf_gauss_7_3 <- gen_hrf(hrf_gaussian, mean = 7, sd = 3, name = "Gaussian (Mean=7, SD=3)")
hrf_gauss_5_2 <- gen_hrf(hrf_gaussian, mean = 5, sd = 2, name = "Gaussian (Mean=5, SD=2)")
hrf_gauss_4_1 <- gen_hrf(hrf_gaussian, mean = 4, sd = 1, name = "Gaussian (Mean=4, SD=1)")

## ----modify_gaussian_params_plot, echo = FALSE--------------------------------
# Compare the HRFs using plot_hrfs()
plot_hrfs(hrf_gauss_7_3, hrf_gauss_5_2, hrf_gauss_4_1,
          title = "Gaussian HRFs with Different Parameters")

## ----blocked_hrfs-------------------------------------------------------------
# Create blocked HRFs using the SPM canonical HRF with different durations
hrf_spm_w1 <- block_hrf(HRF_SPMG1, width = 1)
hrf_spm_w2 <- block_hrf(HRF_SPMG1, width = 2)
hrf_spm_w4 <- block_hrf(HRF_SPMG1, width = 4)

## ----blocked_hrfs_plot, echo = FALSE------------------------------------------
plot_hrfs(hrf_spm_w1, hrf_spm_w2, hrf_spm_w4,
          labels = c("Width=1s", "Width=2s", "Width=4s"),
          title = "SPM Canonical HRF for Different Event Durations",
          subtitle = "Using block_hrf()")

## ----blocked_normalized-------------------------------------------------------
# Create normalized blocked HRFs
hrf_spm_w1_norm <- block_hrf(HRF_SPMG1, width = 1, normalize = TRUE)
hrf_spm_w2_norm <- block_hrf(HRF_SPMG1, width = 2, normalize = TRUE)
hrf_spm_w4_norm <- block_hrf(HRF_SPMG1, width = 4, normalize = TRUE)

## ----blocked_normalized_plot, echo = FALSE------------------------------------
plot_hrfs(hrf_spm_w1_norm, hrf_spm_w2_norm, hrf_spm_w4_norm,
          labels = c("Width=1s", "Width=2s", "Width=4s"),
          title = "Normalized SPM Canonical HRF for Different Durations",
          subtitle = "Using block_hrf(normalize = TRUE)")

## ----blocked_summate_false----------------------------------------------------
# Create non-summating blocked HRFs
hrf_spm_w2_nosum <- block_hrf(HRF_SPMG1, width = 2, summate = FALSE)
hrf_spm_w4_nosum <- block_hrf(HRF_SPMG1, width = 4, summate = FALSE)
hrf_spm_w8_nosum <- block_hrf(HRF_SPMG1, width = 8, summate = FALSE)

## ----blocked_summate_false_plot, echo = FALSE---------------------------------
plot_hrfs(hrf_spm_w2_nosum, hrf_spm_w4_nosum, hrf_spm_w8_nosum,
          labels = c("Width=2s", "Width=4s", "Width=8s"),
          title = "Non-Summating (Saturating) SPM HRF for Different Durations",
          subtitle = "Using block_hrf(summate = FALSE)")

## ----blocked_summate_false_norm-----------------------------------------------
# Create normalized, non-summating blocked HRFs
hrf_spm_w2_nosum_norm <- block_hrf(HRF_SPMG1, width = 2, summate = FALSE, normalize = TRUE)
hrf_spm_w4_nosum_norm <- block_hrf(HRF_SPMG1, width = 4, summate = FALSE, normalize = TRUE)
hrf_spm_w8_nosum_norm <- block_hrf(HRF_SPMG1, width = 8, summate = FALSE, normalize = TRUE)

## ----blocked_summate_false_norm_plot, echo = FALSE----------------------------
plot_hrfs(hrf_spm_w2_nosum_norm, hrf_spm_w4_nosum_norm, hrf_spm_w8_nosum_norm,
          labels = c("Width=2s", "Width=4s", "Width=8s"),
          title = "Normalized, Non-Summating SPM HRF for Different Durations",
          subtitle = "Using block_hrf(summate = FALSE, normalize = TRUE)")

## ----lagged_hrfs--------------------------------------------------------------
# Create lagged versions of the Gaussian HRF
hrf_gauss_lag_neg2 <- lag_hrf(HRF_GAUSSIAN, lag = -2)
hrf_gauss_lag_0 <- HRF_GAUSSIAN # Original (lag=0)
hrf_gauss_lag_pos3 <- lag_hrf(HRF_GAUSSIAN, lag = 3)

## ----lagged_hrfs_plot, echo = FALSE-------------------------------------------
plot_hrfs(hrf_gauss_lag_neg2, hrf_gauss_lag_0, hrf_gauss_lag_pos3,
          labels = c("Lag=-2s", "Lag=0s", "Lag=+3s"),
          title = "Gaussian HRF with Different Temporal Lags",
          subtitle = "Using lag_hrf()")

## ----lagged_blocked_hrfs------------------------------------------------------
# Create HRFs that are both lagged and blocked
hrf_lb_1 <- HRF_GAUSSIAN %>% lag_hrf(1) %>% block_hrf(width = 1, normalize = TRUE)
hrf_lb_3 <- HRF_GAUSSIAN %>% lag_hrf(3) %>% block_hrf(width = 3, normalize = TRUE)
hrf_lb_5 <- HRF_GAUSSIAN %>% lag_hrf(5) %>% block_hrf(width = 5, normalize = TRUE)

## ----lagged_blocked_hrfs_plot, echo = FALSE-----------------------------------
plot_hrfs(hrf_lb_1, hrf_lb_3, hrf_lb_5,
          labels = c("Lag=1, Width=1", "Lag=3, Width=3", "Lag=5, Width=5"),
          title = "Gaussian HRFs with Combined Lag and Duration",
          subtitle = "Using lag_hrf() %>% block_hrf()")

## ----gen_hrf_lag_width--------------------------------------------------------
# Using gen_hrf directly
hrf_lb_gen_3 <- gen_hrf(hrf_gaussian, lag = 3, width = 3, normalize = TRUE)
resp_lb_gen_3 <- hrf_lb_gen_3(time_points)

# Compare (should be very similar to hrf_lb_3 from piped version)
# plot(time_points, resp_lb_3, type = 'l', col = 2, lwd = 2, main = "Piped vs gen_hrf")
# lines(time_points, resp_lb_gen_3, col = 1, lty = 2, lwd = 2)
# legend("topright", legend = c("Piped", "gen_hrf"), col = c(2, 1), lty = c(1, 2), lwd = 2)

## ----spm_basis_sets-----------------------------------------------------------
# SPM + Temporal Derivative (2 basis functions)
print(HRF_SPMG2)

# SPM + Temporal + Dispersion Derivatives (3 basis functions)
print(HRF_SPMG3)

## ----spm_basis_sets_plot, echo = FALSE, fig.height = 4------------------------
resp_spmg2 <- HRF_SPMG2(time_points)
resp_spmg3 <- HRF_SPMG3(time_points)

# Plot SPMG2
matplot(time_points, resp_spmg2, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "SPM + Temporal Derivative Basis Set (HRF_SPMG2)")
legend("topright", legend = c("Canonical", "Temporal Deriv."), col = 1:2, lty = 1, lwd = 1.5)

## ----spm_basis_sets_plot2, echo = FALSE, fig.height = 4-----------------------
# Plot SPMG3
matplot(time_points, resp_spmg3, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "SPM + Temporal + Dispersion Derivative Basis Set (HRF_SPMG3)")
legend("topright", legend = c("Canonical", "Temporal Deriv.", "Dispersion Deriv."), col = 1:3, lty = 1, lwd = 1.5)

## ----bspline_basis------------------------------------------------------------
# B-spline basis with N=5 basis functions, degree=3 (cubic)
hrf_bs_5_3 <- gen_hrf(hrf_bspline, N = 5, degree = 3, name = "B-spline (N=5, deg=3)")
print(hrf_bs_5_3)

# B-spline basis with N=10 basis functions, degree=1 (linear -> tent functions)
hrf_bs_10_1 <- gen_hrf(hrf_bspline, N = 10, degree = 1, name = "Tent Set (N=10)")
print(hrf_bs_10_1)

## ----bspline_basis_plot, echo = FALSE, fig.height = 4-------------------------
resp_bs_5_3 <- hrf_bs_5_3(time_points)
matplot(time_points, resp_bs_5_3, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "B-spline Basis Set (N=5, degree=3)")

## ----tent_basis_plot, echo = FALSE, fig.height = 4----------------------------
resp_bs_10_1 <- hrf_bs_10_1(time_points)
matplot(time_points, resp_bs_10_1, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "Tent Function Basis Set (B-spline, N=10, degree=1)")

## ----sine_basis---------------------------------------------------------------
hrf_sin_5 <- gen_hrf(hrf_sine, N = 5, name = "Sine Basis (N=5)")
print(hrf_sin_5)

## ----sine_basis_plot, echo = FALSE, fig.height = 4----------------------------
resp_sin_5 <- hrf_sin_5(time_points)
matplot(time_points, resp_sin_5, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "Sine Basis Set (N=5)")

## ----half_cosine, echo = FALSE, fig.height = 4--------------------------------
# Use default parameters from Woolrich et al. (2004)
resp_half_cos <- hrf_half_cosine(time_points)
plot(time_points, resp_half_cos, type = 'l', lwd = 1.5,
     xlab = "Time (seconds)", ylab = "BOLD Response",
     main = "Half-Cosine HRF Shape (Woolrich et al., 2004)")

## ----daguerre_basis, eval=FALSE-----------------------------------------------
# # Future implementation:
# # hrf_dag_3 <- HRF_DAGUERRE_BASIS(n_basis = 3)
# # print(hrf_dag_3)
# # resp_dag_3 <- hrf_dag_3(time_points)
# # matplot(time_points, resp_dag_3, type = 'l', lty = 1, lwd = 1.5,
# #         xlab = "Time (seconds)", ylab = "BOLD Response",
# #         main = "Daguerre Basis Set (N=3)")
# ``` -->
# 
# ## Other HRF Shapes
# 
# ### Gamma HRF
# 
# The `hrf_gamma` function uses the gamma probability density function.
# 

## ----gamma_hrf----------------------------------------------------------------
hrf_gam <- gen_hrf(hrf_gamma, shape = 6, rate = 1, name = "Gamma (shape=6, rate=1)")
print(hrf_gam)

## ----gamma_hrf_plot, echo = FALSE, fig.height = 4-----------------------------
plot(hrf_gam, time = time_points)

## ----mexhat_hrf---------------------------------------------------------------
hrf_mh <- gen_hrf(hrf_mexhat, mean = 6, sd = 1.5, name = "Mexican Hat (mean=6, sd=1.5)")
print(hrf_mh)

## ----mexhat_hrf_plot, echo = FALSE, fig.height = 4----------------------------
plot(hrf_mh, time = time_points)

## ----inv_logit_hrf------------------------------------------------------------
hrf_il <- gen_hrf(hrf_inv_logit, mu1 = 5, s1 = 1, mu2 = 15, s2 = 1.5, name = "Inv. Logit Diff.")
print(hrf_il)

## ----inv_logit_hrf_plot, echo = FALSE, fig.height = 4-------------------------
plot(hrf_il, time = time_points)

## ----boxcar_basic-------------------------------------------------------------
# Create a boxcar of width 5 seconds (from 0 to 5 seconds)
hrf_box <- hrf_boxcar(width = 5)
print(hrf_box)

## ----boxcar_basic_plot, echo = FALSE, fig.height = 4--------------------------
# Evaluate
t <- seq(-2, 10, by = 0.1)
resp_box <- evaluate(hrf_box, t)

plot(t, resp_box, type = 's', lwd = 2,
     xlab = "Time (seconds)", ylab = "Response",
     main = "Simple Boxcar HRF (width = 5 seconds)")
abline(h = 0, lty = 2, col = "gray")

## ----boxcar_delayed-----------------------------------------------------------
# Boxcar from 4-8 seconds post-stimulus (capturing the expected BOLD peak)
# Use lag_hrf() to delay a 4-second boxcar by 4 seconds
hrf_delayed <- hrf_boxcar(width = 4) %>% lag_hrf(lag = 4)

## ----boxcar_delayed_plot, echo = FALSE----------------------------------------
# Compare boxcar with normalized SPM canonical
plot_hrfs(hrf_delayed, HRF_SPMG1,
          labels = c("Boxcar (4-8s)", "SPM Canonical"),
          normalize = TRUE,
          title = "Boxcar vs. Traditional HRF",
          subtitle = "Boxcar captures a specific time window; traditional HRF models hemodynamic delay")

## ----boxcar_normalized--------------------------------------------------------
# Normalized boxcar - integral = 1
# A 4-second boxcar lagged by 4 seconds (captures 4-8s window)
hrf_norm <- hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 4)

# Check: amplitude should be 1/4 = 0.25
t_fine <- seq(0, 12, by = 0.01)
resp_norm <- evaluate(hrf_norm, t_fine)
cat("Amplitude of normalized boxcar:", max(resp_norm), "\n")
cat("Expected (1/width):", 1/4, "\n")

# Verify integral ≈ 1
integral <- sum(resp_norm) * 0.01
cat("Integral of normalized boxcar:", round(integral, 3), "\n")

## ----weighted_width-----------------------------------------------------------
# 6 weights evenly spaced over 10 seconds (at 0, 2, 4, 6, 8, 10)
hrf_wt_width <- hrf_weighted(
  weights = c(0.1, 0.3, 1.0, 1.0, 0.3, 0.1),
  width = 10,
  method = "constant"
)

## ----weighted_width_plot, echo = FALSE, fig.height = 4------------------------
resp_wt_width <- evaluate(hrf_wt_width, time_points)

plot(time_points, resp_wt_width, type = 's', lwd = 2,
     xlab = "Time (seconds)", ylab = "Weight",
     main = "Weighted Step Function HRF (width = 10)")
abline(h = 0, lty = 2, col = "gray")

## ----weighted_times-----------------------------------------------------------
# Weighted step function with explicit time points
hrf_wt <- hrf_weighted(
  weights = c(0.1, 0.3, 1.0, 1.0, 0.3, 0.1),
  times = c(2, 4, 6, 8, 10, 12),
  method = "constant"
)

## ----weighted_times_plot, echo = FALSE, fig.height = 4------------------------
resp_wt <- evaluate(hrf_wt, time_points)

plot(time_points, resp_wt, type = 's', lwd = 2,
     xlab = "Time (seconds)", ylab = "Weight",
     main = "Weighted Step Function HRF (explicit times)")
abline(h = 0, lty = 2, col = "gray")

## ----weighted_linear----------------------------------------------------------
# Smooth weights using linear interpolation
hrf_smooth <- hrf_weighted(
  weights = c(0, 0.3, 1.0, 1.0, 0.3, 0),
  times = c(2, 4, 6, 8, 10, 12),
  method = "linear"
)

## ----weighted_linear_plot, echo = FALSE---------------------------------------
resp_smooth <- evaluate(hrf_smooth, time_points)

# Compare constant vs. linear interpolation
plot_df_wt <- data.frame(
  Time = time_points,
  `Step (constant)` = resp_wt,
  `Smooth (linear)` = resp_smooth
) %>%
  pivot_longer(-Time, names_to = "Method", values_to = "Response")

ggplot(plot_df_wt, aes(x = Time, y = Response, color = Method)) +
  geom_line(linewidth = 1) +
  labs(title = "Weighted HRF: Constant vs. Linear Interpolation",
       x = "Time (seconds)",
       y = "Weight") +
  theme_minimal()

## ----weighted_subsecond-------------------------------------------------------
# Sub-second intervals: create a Gaussian-shaped weight function
times_fine <- seq(4, 10, by = 0.25)
weights_gaussian <- dnorm(times_fine, mean = 7, sd = 1)

hrf_gauss_wt <- hrf_weighted(weights_gaussian, times = times_fine, method = "linear")

## ----weighted_subsecond_plot, echo = FALSE, fig.height = 4--------------------
resp_gauss_wt <- evaluate(hrf_gauss_wt, time_points)

plot(time_points, resp_gauss_wt, type = 'l', lwd = 2,
     xlab = "Time (seconds)", ylab = "Weight",
     main = "Gaussian-Shaped Weights (Sub-second Precision)")

## ----weighted_normalized------------------------------------------------------
# Normalized weights - creates weighted average interpretation
hrf_wt_norm <- hrf_weighted(
  weights = c(1, 2, 2, 1),  # Will be normalized
  times = c(4, 6, 8, 10),
  method = "constant",
  normalize = TRUE
)

# The coefficient β will estimate: (1*Y[4-6] + 2*Y[6-8] + 2*Y[8-10] + 1*Y[10+]) / 6
# where Y[a-b] is the signal in that interval

t_check <- seq(0, 12, by = 0.01)
resp_wt_norm <- evaluate(hrf_wt_norm, t_check)

# Verify: integral should be approximately 1
integral_wt <- sum(resp_wt_norm) * 0.01
cat("Integral of normalized weighted HRF:", round(integral_wt, 3), "\n")

## ----early_late_comparison----------------------------------------------------
# Early window: 2-6 seconds (4-second boxcar lagged by 2 seconds)
hrf_early <- hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 2)

# Late window: 8-12 seconds (4-second boxcar lagged by 8 seconds)
hrf_late <- hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 8)

## ----early_late_comparison_plot, echo = FALSE---------------------------------
# Evaluate both
resp_early <- evaluate(hrf_early, time_points)
resp_late <- evaluate(hrf_late, time_points)

# Also show the canonical HRF for reference
resp_spm_ref <- HRF_SPMG1(time_points)
resp_spm_ref <- resp_spm_ref / max(resp_spm_ref) * 0.3  # Scale for visibility

plot_df_windows <- data.frame(
  Time = time_points,
  `Early (2-6s)` = resp_early,
  `Late (8-12s)` = resp_late,
  `SPM (scaled)` = resp_spm_ref
) %>%
  pivot_longer(-Time, names_to = "Window", values_to = "Response")

ggplot(plot_df_windows, aes(x = Time, y = Response, color = Window)) +
  geom_line(linewidth = 1) +
  labs(title = "Early vs. Late Response Windows",
       subtitle = "Normalized boxcars for extracting mean signal in different windows",
       x = "Time (seconds)",
       y = "Response") +
  theme_minimal() +
  scale_color_manual(values = c("Early (2-6s)" = "blue",
                                 "Late (8-12s)" = "red",
                                 "SPM (scaled)" = "gray50"))

## ----boxcar_regressor---------------------------------------------------------
# Create a regressor with boxcar HRF (4-second window starting 4s after onset)
reg_boxcar <- regressor(
  onsets = c(0, 20, 40),
  hrf = hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 4)
)

# Compare with traditional SPM HRF
reg_spm <- regressor(onsets = c(0, 20, 40), hrf = HRF_SPMG1)

## ----boxcar_regressor_plot, echo = FALSE--------------------------------------
# Evaluate the design regressor
t_design <- seq(0, 60, by = 0.5)
design_boxcar <- evaluate(reg_boxcar, t_design)

design_spm <- evaluate(reg_spm, t_design)
design_spm <- design_spm / max(design_spm)  # Normalize for comparison

plot_df_design <- data.frame(
  Time = t_design,
  `Boxcar (4-8s)` = design_boxcar,
  `SPM Canonical` = design_spm
) %>%
  pivot_longer(-Time, names_to = "HRF", values_to = "Response")

ggplot(plot_df_design, aes(x = Time, y = Response, color = HRF)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Design Matrix Regressors: Boxcar vs. Traditional HRF",
       subtitle = "Events at t = 0, 20, 40 seconds",
       x = "Time (seconds)",
       y = "Regressor Value") +
  theme_minimal()

## ----custom_basis_lagged------------------------------------------------------
# Create a list of lagged Gaussian HRFs
lag_times <- seq(0, 10, by = 2)
list_of_hrfs <- lapply(lag_times, function(lag) {
  lag_hrf(HRF_GAUSSIAN, lag = lag)
})

# Combine them into a single HRF basis set object
hrf_custom_set <- do.call(gen_hrf_set, list_of_hrfs)
print(hrf_custom_set) # Note: name is default 'hrf_set', nbasis is 6

## ----custom_basis_lagged_plot, echo = FALSE, fig.height = 4-------------------
# Evaluate and plot
resp_custom_set <- hrf_custom_set(time_points)
matplot(time_points, resp_custom_set, type = 'l', lty = 1, lwd = 1.5,
        xlab = "Time (seconds)", ylab = "BOLD Response",
        main = "Custom Basis Set (Lagged Gaussians)")

## ----empirical_hrf_single-----------------------------------------------------
# Simulate an average measured response profile
sim_times <- 0:24
set.seed(42) # For reproducibility
sim_profile <- rowMeans(replicate(20, {
  h <- HRF_SPMG1 %>% lag_hrf(lag = runif(n = 1, min = -1, max = 1)) %>%
                    block_hrf(width = runif(n = 1, min = 0, max = 2))
  h(sim_times)
}))

# Normalize profile to max = 1 for better visualization
sim_profile_norm <- sim_profile / max(sim_profile)

# Create the empirical HRF function from the normalized profile
emp_hrf <- gen_empirical_hrf(sim_times, sim_profile_norm)
print(emp_hrf)

## ----empirical_hrf_single_plot, echo = FALSE----------------------------------
# Evaluate and plot (using a finer time grid for interpolation)
fine_times <- seq(0, 24, by = 0.1)
resp_emp <- emp_hrf(fine_times)

# Plot the interpolated curve with the original points
plot(fine_times, resp_emp, type = 'l', lwd = 1.5,
     xlab = "Time (seconds)", ylab = "BOLD Response",
     main = "Empirical HRF from Simulated Average Profile")
points(sim_times, sim_profile_norm, pch = 16, col = "red", cex = 1) # Show original points

## ----empirical_hrf_pca--------------------------------------------------------
# 1. Simulate a matrix of diverse HRFs
set.seed(123) # for reproducibility
n_sim <- 50
sim_mat <- replicate(n_sim, {
  hrf_func <- HRF_SPMG1 %>%
              lag_hrf(lag = runif(1, -2, 2)) %>%
              block_hrf(width = runif(1, 0, 3))
  hrf_func(sim_times)
})

## ----empirical_hrf_pca_plot1, echo = FALSE------------------------------------
# Show a sample of simulated HRFs to illustrate variability
matplot(sim_times, sim_mat[, 1:10], type = 'l', col = scales::alpha("gray", 0.7), lty = 1,
        xlab = "Time (seconds)", ylab = "Response",
        main = "Sample of Simulated HRF Profiles")

## ----empirical_hrf_pca2-------------------------------------------------------
# 2. Perform PCA on the transpose (each column = one HRF, each row = one time point)
pca_res <- prcomp(t(sim_mat), center = TRUE, scale. = FALSE)
n_components <- 3

# Print variance explained by top components
variance_explained <- summary(pca_res)$importance[2, 1:n_components]
cat("Variance explained by top", n_components, "components:",
    paste0(round(variance_explained * 100, 1), "%"), "\n")

# Extract the top principal components
pc_vectors <- pca_res$rotation[, 1:n_components]

# 3. Convert principal components into HRF functions
list_pc_hrfs <- list()

for (i in 1:n_components) {
  pc_vec <- pc_vectors[, i]
  pc_vec_zeroed <- pc_vec - pc_vec[1]
  max_abs <- max(abs(pc_vec_zeroed))
  pc_vec_norm <- pc_vec_zeroed / max_abs
  list_pc_hrfs[[i]] <- gen_empirical_hrf(sim_times, pc_vec_norm)
}

# 4. Combine PC HRFs into a basis set using gen_hrf_set
emp_pca_basis <- do.call(gen_hrf_set, list_pc_hrfs)
print(emp_pca_basis)

## ----empirical_hrf_pca_plot2, echo = FALSE------------------------------------
# Evaluate and plot the basis functions
resp_pca_basis <- emp_pca_basis(sim_times)

pc_df <- as.data.frame(resp_pca_basis)
names(pc_df) <- paste("PC", 1:n_components)
pc_df$Time <- sim_times

pc_df_long <- pivot_longer(pc_df, -Time, names_to = "Component", values_to = "Value")

ggplot(pc_df_long, aes(x = Time, y = Value, color = Component)) +
  geom_line(linewidth = 1.2) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Empirical Basis Set from PCA",
       subtitle = paste0("First ", n_components, " Principal Components"),
       x = "Time (seconds)",
       y = "Component Value") +
  theme_minimal() +
  theme(legend.position = "right")

