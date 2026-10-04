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
library(dplyr) # For pipe operator %>%
library(ggplot2) # For plotting
library(tidyr) # For data manipulation

## ----print_hrfs--------------------------------------------------
# SPM canonical HRF (based on difference of two gamma functions)
print(HRF_SPMG1)

# Gaussian HRF
print(HRF_GAUSSIAN)

## ----evaluate_basic_hrfs, fig.alt="SPM and Gaussian HRFs normalized to a peak of one. The SPM response includes a negative undershoot after about 12 seconds; the Gaussian does not."----
time_points <- seq(0, 25, by = 0.1)

# normalize = TRUE scales each HRF to peak at 1.0
plot_hrfs(HRF_SPMG1, HRF_GAUSSIAN,
          labels = c("SPM canonical", "Gaussian"),
          normalize = TRUE,
          time = time_points,
          title = "SPM and Gaussian HRFs",
          subtitle = "Each curve normalized to peak at 1")

## ----fixed-hrf-scale---------------------------------------------
spm_scaled <- gen_hrf(HRF_SPMG1, hrf_norm = "spm")
spm_grid <- seq(0, 32, length.out = 1600)
sum(spm_scaled(spm_grid))

## ----check-fixed-hrf-scale, include = FALSE----------------------
stopifnot(abs(sum(spm_scaled(spm_grid)) - 1) < 1e-12)

## ----modify_gaussian_params--------------------------------------
# Create Gaussian HRFs with different parameters using gen_hrf
# Note: hrf_gaussian is the underlying function, not the HRF object HRF_GAUSSIAN
hrf_gauss_4_1 <- gen_hrf(hrf_gaussian, mean = 4, sd = 1, name = "Gaussian (Mean=4, SD=1)")
hrf_gauss_5_2 <- gen_hrf(hrf_gaussian, mean = 5, sd = 2, name = "Gaussian (Mean=5, SD=2)")
hrf_gauss_7_3 <- gen_hrf(hrf_gaussian, mean = 7, sd = 3, name = "Gaussian (Mean=7, SD=3)")

## ----modify_gaussian_params_plot, fig.alt="Three Gaussian HRFs. Later means peak later, and larger standard deviations give lower, broader curves."----
plot_hrfs(hrf_gauss_4_1, hrf_gauss_5_2, hrf_gauss_7_3,
          labels = c("mean 4, sd 1", "mean 5, sd 2", "mean 7, sd 3"),
          time = time_points, palette = "ordered",
          title = "Gaussian HRF parameters",
          subtitle = "Later mean, later peak; larger sd, broader")

## ----blocked_hrfs------------------------------------------------
# Create blocked HRFs using the SPM canonical HRF with different durations
hrf_spm_w1 <- block_hrf(HRF_SPMG1, width = 1)
hrf_spm_w2 <- block_hrf(HRF_SPMG1, width = 2)
hrf_spm_w4 <- block_hrf(HRF_SPMG1, width = 4)

## ----blocked_hrfs_plot, fig.alt="SPM canonical responses to 1, 2 and 4 second events. Longer events produce higher and later peaks."----
plot_hrfs(hrf_spm_w1, hrf_spm_w2, hrf_spm_w4,
          labels = c("1 s event", "2 s event", "4 s event"),
          time = time_points, palette = "ordered",
          title = "Longer events: larger responses",
          subtitle = "block_hrf(HRF_SPMG1, width = 1, 2, 4)")

## ----blocked_normalized------------------------------------------
# Create normalized blocked HRFs
hrf_spm_w1_norm <- block_hrf(HRF_SPMG1, width = 1, normalize = TRUE)
hrf_spm_w2_norm <- block_hrf(HRF_SPMG1, width = 2, normalize = TRUE)
hrf_spm_w4_norm <- block_hrf(HRF_SPMG1, width = 4, normalize = TRUE)

## ----blocked_normalized_plot, fig.alt="Normalized SPM responses to 1, 2 and 4 second events. All peak at 1; longer events peak later and are broader."----
plot_hrfs(hrf_spm_w1_norm, hrf_spm_w2_norm, hrf_spm_w4_norm,
          labels = c("1 s event", "2 s event", "4 s event"),
          time = time_points, palette = "ordered",
          title = "Normalized: only shape changes",
          subtitle = "block_hrf(..., normalize = TRUE)")

## ----blocked_summate_false---------------------------------------
# Create non-summating blocked HRFs
hrf_spm_w2_nosum <- block_hrf(HRF_SPMG1, width = 2, summate = FALSE)
hrf_spm_w4_nosum <- block_hrf(HRF_SPMG1, width = 4, summate = FALSE)
hrf_spm_w8_nosum <- block_hrf(HRF_SPMG1, width = 8, summate = FALSE)

## ----blocked_summate_false_plot, fig.alt="Averaged (summate = FALSE) SPM responses to 2, 4 and 8 second events. Peaks do not grow with duration; the 8 second response is lower and broader."----
plot_hrfs(hrf_spm_w2_nosum, hrf_spm_w4_nosum, hrf_spm_w8_nosum,
          labels = c("2 s event", "4 s event", "8 s event"),
          time = time_points, palette = "ordered",
          title = "Averaging: peaks do not grow",
          subtitle = "block_hrf(..., summate = FALSE)")

## ----lagged_hrfs-------------------------------------------------
# Create lagged versions of the Gaussian HRF
hrf_gauss_lag_neg2 <- lag_hrf(HRF_GAUSSIAN, lag = -2)
hrf_gauss_lag_0 <- HRF_GAUSSIAN # Original (lag=0)
hrf_gauss_lag_pos3 <- lag_hrf(HRF_GAUSSIAN, lag = 3)

## ----lagged_hrfs_plot, fig.alt="Gaussian HRF shifted by -2, 0 and +3 seconds. The shape is unchanged; only the peak time moves."----
plot_hrfs(hrf_gauss_lag_neg2, hrf_gauss_lag_0, hrf_gauss_lag_pos3,
          labels = c("lag -2 s", "lag 0 s", "lag +3 s"),
          time = time_points, palette = "ordered",
          title = "Lag moves the response in time",
          subtitle = "lag_hrf(HRF_GAUSSIAN, lag = -2, 0, 3)")

## ----lagged_blocked_hrfs-----------------------------------------
# Create HRFs that are both lagged and blocked
hrf_lb_1 <- HRF_GAUSSIAN %>% lag_hrf(1) %>% block_hrf(width = 1, normalize = TRUE)
hrf_lb_3 <- HRF_GAUSSIAN %>% lag_hrf(3) %>% block_hrf(width = 3, normalize = TRUE)
hrf_lb_5 <- HRF_GAUSSIAN %>% lag_hrf(5) %>% block_hrf(width = 5, normalize = TRUE)

## ----lagged_blocked_hrfs_plot, fig.alt="Gaussian HRFs with lag and width both equal to 1, 3 and 5 seconds, normalized to peak at 1. Larger values peak later and are broader."----
plot_hrfs(hrf_lb_1, hrf_lb_3, hrf_lb_5,
          labels = c("lag 1 s, width 1 s", "lag 3 s, width 3 s", "lag 5 s, width 5 s"),
          time = time_points, palette = "ordered",
          title = "Lag and duration combined",
          subtitle = "lag_hrf() %>% block_hrf(normalize = TRUE)")

## ----gen_hrf_lag_width-------------------------------------------
hrf_lb_gen_3 <- gen_hrf(hrf_gaussian, lag = 3, width = 3, normalize = TRUE)

# Largest difference from the piped version over 0-25 s
max(abs(hrf_lb_gen_3(time_points) - hrf_lb_3(time_points)))

## ----spm_basis_sets----------------------------------------------
# SPM + Temporal Derivative (2 basis functions)
print(HRF_SPMG2)

# SPM + Temporal + Dispersion Derivatives (3 basis functions)
print(HRF_SPMG3)

## ----spm_basis_sets_plot, fig.height = 4.4, fig.alt="Canonical SPM response and its temporal and dispersion derivatives, with all positive and negative values retained."----
plot_hrfs(HRF_SPMG3,
          labels = c("Canonical", "Temporal derivative", "Dispersion derivative"),
          time = time_points,
          title = "SPM canonical and derivatives")

## ----spm_derivative_shift, fig.height = 4.4, fig.alt="Canonical SPM HRF (dashed grey) and the canonical plus or minus 0.5 times its temporal derivative. Adding the derivative moves the peak earlier; subtracting it moves the peak later."----
spm_early <- gen_empirical_hrf(time_points,
  drop(HRF_SPMG2(time_points) %*% c(1, 0.5)))
spm_late <- gen_empirical_hrf(time_points,
  drop(HRF_SPMG2(time_points) %*% c(1, -0.5)))

plot_hrfs(spm_early, spm_late,
          labels = c("canonical + 0.5 x derivative", "canonical - 0.5 x derivative"),
          time = time_points, reference = HRF_SPMG1, reference_label = "canonical",
          title = "The derivative shifts the peak")

## ----bspline_basis-----------------------------------------------
# B-spline basis with N=5 basis functions, degree=3 (cubic)
hrf_bs_5_3 <- gen_hrf(hrf_bspline, N = 5, degree = 3, name = "B-spline (N=5, deg=3)")
print(hrf_bs_5_3)

# B-spline basis with N=11 basis functions, degree=1 (linear -> tent functions):
# tents peak every 2 s, at 2, 4, ..., 22 s
hrf_bs_10_1 <- gen_hrf(hrf_bspline, N = 11, degree = 1, name = "Tent Set (N=11)")
print(hrf_bs_10_1)

## ----bspline_basis_plot, fig.alt="Five cubic B-spline basis functions on 0 to 24 seconds, labelled B1 to B5 at their peaks, tiling the window from early to late."----
bspline_times <- seq(0, 24, by = 0.1)
plot_hrfs(hrf_bs_5_3, time = bspline_times,
          title = "Cubic B-spline basis (N = 5)")

## ----tent_basis_plot, fig.alt="Eleven piecewise-linear tent functions peaking every 2 seconds from 2 to 22 seconds, each overlapping its neighbours and all zero at 0 and 24 seconds."----
plot_hrfs(hrf_bs_10_1, time = bspline_times,
          title = "Tent basis: linear B-splines (N = 11)")

## ----sine_basis--------------------------------------------------
hrf_sin_5 <- gen_hrf(hrf_sine, N = 5, name = "Sine Basis (N=5)")
print(hrf_sin_5)

## ----sine_basis_plot, fig.height = 5, fig.alt="Five sine basis functions in separate panels, with one to five full cycles over the 24 second window."----
plot_hrfs(hrf_sin_5, time = bspline_times, layout = "stack",
          title = "Sine basis (N = 5)")

## ----half_cosine, fig.alt="Two half-cosine HRFs. The default has no dip or undershoot; with f1 = -0.1 and f2 = -0.2 an initial dip to -0.1 and an undershoot to -0.2 at about 13 seconds appear."----
hrf_hc_default <- gen_hrf(hrf_half_cosine, name = "Half-cosine (default)")
hrf_hc_dip <- gen_hrf(hrf_half_cosine, f1 = -0.1, f2 = -0.2,
                      name = "Half-cosine (f1 = -0.1, f2 = -0.2)")
plot_hrfs(hrf_hc_default, hrf_hc_dip,
          labels = c("defaults (f1 = f2 = 0)", "f1 = -0.1, f2 = -0.2"),
          time = time_points,
          title = "Half-cosine HRF",
          subtitle = "Woolrich et al. (2004)")

## ----other_shapes------------------------------------------------
# Gamma probability density
hrf_gam <- gen_hrf(hrf_gamma, shape = 6, rate = 1, name = "Gamma (shape=6, rate=1)")

# Mexican hat wavelet (second derivative of a Gaussian)
hrf_mh <- gen_hrf(hrf_mexhat, mean = 6, sd = 1.5, name = "Mexican Hat (mean=6, sd=1.5)")

# Difference of two inverse-logit (sigmoid) functions: separate rise and fall
hrf_il <- gen_hrf(hrf_inv_logit, mu1 = 5, s1 = 1, mu2 = 15, s2 = 1.5, name = "Inv. Logit Diff.")

## ----other_shapes_plot, fig.height = 5, fig.alt="Three HRF shapes in separate panels: a gamma peaking at 5 seconds, a Mexican hat with negative side lobes around a peak at 6 seconds, and an inverse-logit difference that rises around 5 seconds and falls around 15 seconds."----
plot_hrfs(hrf_gam, hrf_mh, hrf_il,
          labels = c("Gamma (shape 6, rate 1)", "Mexican hat (mean 6, sd 1.5)",
                     "Inverse-logit (rise 5 s, fall 15 s)"),
          time = time_points, layout = "stack",
          title = "Other HRF shapes")

## ----boxcar_basic------------------------------------------------
# Create a boxcar of width 5 seconds (from 0 to 5 seconds)
hrf_box <- hrf_boxcar(width = 5)
print(hrf_box)

## ----boxcar_delayed----------------------------------------------
# Boxcar from 4-8 seconds post-stimulus (capturing the expected BOLD peak)
# Use lag_hrf() to delay a 4-second boxcar by 4 seconds
hrf_delayed <- hrf_boxcar(width = 4) %>% lag_hrf(lag = 4)

## ----boxcar_delayed_plot, fig.alt="A 0 to 5 second boxcar, a 4 to 8 second boxcar, and the normalized SPM canonical HRF as a dashed grey reference. The delayed boxcar covers the rising edge and peak of the canonical response."----
plot_hrfs(hrf_box, hrf_delayed,
          labels = c("Boxcar, 0-5 s", "Boxcar, 4-8 s"),
          normalize = TRUE, time = seq(0, 25, by = 0.05),
          reference = HRF_SPMG1, reference_label = "SPM canonical (peak = 1)",
          title = "Boxcars select a time window",
          subtitle = "No hemodynamic delay: a boxcar weights only its window")

## ----boxcar_normalized-------------------------------------------
# Unit-area boxcar: a 4-second window lagged by 4 seconds (4-8 s)
hrf_norm <- hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 4)

# Check: amplitude should be 1/4 = 0.25
t_fine <- seq(0, 12, by = 0.01)
resp_norm <- evaluate(hrf_norm, t_fine)
cat("Amplitude of normalized boxcar:", max(resp_norm), "\n")
cat("Expected (1/width):", 1/4, "\n")

# Verify integral ≈ 1
integral <- sum(resp_norm) * 0.01
cat("Integral of normalized boxcar:", round(integral, 3), "\n")

## ----boxcar_beta_check-------------------------------------------
scan_t <- seq(0, 40, by = 0.5)
y <- ifelse(scan_t >= 4 & scan_t < 8, 5, 0)
beta <- function(hrf) {
  x <- evaluate(regressor(0, hrf), scan_t)
  unname(coef(lm(y ~ x))["x"])
}
c(height_1 = beta(hrf_boxcar(width = 4) %>% lag_hrf(lag = 4)),  # about 5: the mean
  unit_area = beta(hrf_norm))                                   # about 20: 5 x 4 s

## ----weighted_width----------------------------------------------
# 6 weights over a 10-second window: six equal bins of 10/6 s
hrf_wt_width <- hrf_weighted(
  weights = c(0.1, 0.3, 1.0, 1.0, 0.3, 0.1),
  width = 10,
  method = "constant"
)

# The same weights in 2-second bins starting at 2, 4, ..., 12 s (window 2-14 s)
hrf_wt <- hrf_weighted(
  weights = c(0.1, 0.3, 1.0, 1.0, 0.3, 0.1),
  times = c(2, 4, 6, 8, 10, 12),
  method = "constant"
)

# Smooth weights using linear interpolation between the same time points
hrf_smooth <- hrf_weighted(
  weights = c(0.1, 0.3, 1.0, 1.0, 0.3, 0.1),
  times = c(2, 4, 6, 8, 10, 12),
  method = "linear"
)

## ----weighted_plot, fig.height = 5.2, fig.alt="Three weighted HRFs in separate panels: steps from 0 to 10 seconds, the same steps from 2 to 12 seconds, and a linear interpolation through the same weights from 2 to 12 seconds."----
plot_hrfs(hrf_wt_width, hrf_wt, hrf_smooth,
          labels = c("width = 10: six bins, 0-10 s",
                     "times = 2, ..., 12: bins, 2-14 s",
                     "method = linear: 2-12 s"),
          time = seq(0, 16, by = 0.02), layout = "stack",
          title = "Weighted HRFs", subtitle = "Weights 0.1, 0.3, 1, 1, 0.3, 0.1")

## ----weighted_subsecond, fig.alt="A Gaussian-shaped weighting function centred at 7 seconds, built from weights every 0.25 seconds between 4 and 10 seconds."----
# Sub-second intervals: create a Gaussian-shaped weight function
times_fine <- seq(4, 10, by = 0.25)
weights_gaussian <- dnorm(times_fine, mean = 7, sd = 1)

hrf_gauss_wt <- hrf_weighted(weights_gaussian, times = times_fine, method = "linear")

plot_hrfs(hrf_gauss_wt, labels = "Gaussian weights, 0.25 s spacing",
          time = seq(0, 14, by = 0.02),
          title = "Weights at sub-second spacing")

## ----weighted_normalized-----------------------------------------
hrf_wt_norm <- hrf_weighted(
  weights = c(1, 2, 2, 1),  # Will be normalized
  times = c(4, 6, 8, 10),
  method = "constant",
  normalize = TRUE
)

# Four 2-second bins (4-6, 6-8, 8-10, 10-12 s); weights 1, 2, 2, 1 are
# rescaled to 1/6, 2/6, 2/6, 1/6
evaluate(hrf_wt_norm, c(5, 7, 9, 11, 13))

## ----early_late_comparison---------------------------------------
# Early window: 2-6 seconds (4-second boxcar lagged by 2 seconds)
hrf_early <- hrf_boxcar(width = 4) %>% lag_hrf(lag = 2)

# Late window: 8-12 seconds (4-second boxcar lagged by 8 seconds)
hrf_late <- hrf_boxcar(width = 4) %>% lag_hrf(lag = 8)

## ----early_late_comparison_plot, fig.alt="Early (2 to 6 seconds) and late (8 to 12 seconds) boxcar windows of height 1, drawn with the SPM canonical HRF scaled to the same height. The early window covers the peak; the late window covers the decline."----
spm_ref <- gen_empirical_hrf(time_points,
  HRF_SPMG1(time_points) / max(HRF_SPMG1(time_points)))
plot_hrfs(hrf_early, hrf_late,
          labels = c("Early window, 2-6 s", "Late window, 8-12 s"),
          time = seq(0, 25, by = 0.05),
          reference = spm_ref, reference_label = "SPM canonical, scaled to peak 1",
          title = "Early and late windows",
          subtitle = "Each beta is the mean signal in its window")

## ----boxcar_regressor--------------------------------------------
# Create a regressor with boxcar HRF (4-second window starting 4s after onset)
reg_boxcar <- regressor(
  onsets = c(0, 20, 40),
  hrf = hrf_boxcar(width = 4, normalize = TRUE) %>% lag_hrf(lag = 4)
)

# Compare with traditional SPM HRF
reg_spm <- regressor(onsets = c(0, 20, 40), hrf = HRF_SPMG1)

## ----boxcar_regressor_plot, fig.height = 4.4, fig.alt="Two regressors for events at 0, 20 and 40 seconds, in separate panels: the smooth SPM canonical regressor, and a boxcar regressor with 4 second plateaus from 4 to 8 seconds after each event."----
plot_regressors(reg_spm, reg_boxcar,
                labels = c("SPM canonical HRF", "Boxcar HRF, 4-8 s window"),
                grid = seq(0, 60, by = 0.05), layout = "stack",
                title = "Boxcar vs. canonical regressor",
                subtitle = "Events at t = 0, 20, 40 s")

## ----custom_basis_lagged-----------------------------------------
# Create a list of lagged Gaussian HRFs
lag_times <- seq(0, 10, by = 2)
list_of_hrfs <- lapply(lag_times, function(lag) {
  lag_hrf(HRF_GAUSSIAN, lag = lag)
})

# Combine them into a single HRF basis set object
hrf_custom_set <- do.call(hrf_set, list_of_hrfs)
print(hrf_custom_set) # Note: name is default 'hrf_set', nbasis is 6

## ----custom_basis_lagged_plot, fig.alt="Six Gaussian basis functions, labelled B1 to B6, peaking every 2 seconds from 6 to 16 seconds."----
plot_hrfs(hrf_custom_set, time = time_points,
          title = "Lagged-Gaussian basis set",
          subtitle = "HRF_GAUSSIAN lagged by 0, 2, ..., 10 s")

## ----empirical_hrf_single----------------------------------------
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

## ----empirical_hrf_single_plot, echo = FALSE, fig.alt="Empirical HRF drawn as a line through 25 measured points, one per second, peaking near 6 seconds with an undershoot after 12 seconds."----
emp_df <- data.frame(time = seq(0, 24, by = 0.1))
emp_df$response <- emp_hrf(emp_df$time)
ggplot(emp_df, aes(time, response)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(colour = hrf_palette(1), linewidth = 0.9) +
  geom_point(data = data.frame(time = sim_times, response = sim_profile_norm),
             colour = hrf_palette(1), size = 1.8) +
  labs(title = "Empirical HRF from a profile",
       x = "Time (s)", y = "Response / peak") +
  theme(plot.title.position = "plot")

## ----empirical_hrf_pca-------------------------------------------
# 1. Simulate a matrix of diverse HRFs
set.seed(123) # for reproducibility
n_sim <- 50
sim_mat <- replicate(n_sim, {
  hrf_func <- HRF_SPMG1 %>%
              lag_hrf(lag = runif(1, -2, 2)) %>%
              block_hrf(width = runif(1, 0, 3))
  hrf_func(sim_times)
})

## ----empirical_hrf_pca_plot1, echo = FALSE, fig.alt="Fifty simulated HRFs drawn as thin grey lines, with their mean as a thick line. Peak times range from about 3 to 8 seconds."----
sim_df <- data.frame(time = rep(sim_times, n_sim),
                     response = as.vector(sim_mat),
                     id = rep(seq_len(n_sim), each = length(sim_times)))
ggplot(sim_df, aes(time, response)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(group = id), colour = "grey60", linewidth = 0.3, alpha = 0.6) +
  geom_line(data = data.frame(time = sim_times, response = rowMeans(sim_mat)),
            colour = hrf_palette(1), linewidth = 1.1) +
  labs(title = "50 simulated HRFs and their mean",
       subtitle = "Random lag and duration",
       x = "Time (s)", y = "Response") +
  theme(plot.title.position = "plot")

## ----empirical_hrf_pca2------------------------------------------
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
  # The sign of a principal component is arbitrary; make its largest
  # deviation positive so PC1 looks like an HRF rather than its mirror image
  pc_vec <- pc_vec * sign(pc_vec[which.max(abs(pc_vec - pc_vec[1]))] - pc_vec[1])
  pc_vec_zeroed <- pc_vec - pc_vec[1]
  max_abs <- max(abs(pc_vec_zeroed))
  pc_vec_norm <- pc_vec_zeroed / max_abs
  list_pc_hrfs[[i]] <- gen_empirical_hrf(sim_times, pc_vec_norm)
}

# 4. Combine PC HRFs into a basis set using hrf_set
emp_pca_basis <- do.call(hrf_set, list_pc_hrfs)
print(emp_pca_basis)

## ----empirical_hrf_pca_plot2, fig.alt="The first three principal components of the simulated HRFs, each scaled to a maximum absolute value of one. PC1 resembles the average HRF; PC2 and PC3 have alternating positive and negative lobes that shift and reshape the response."----
plot_hrfs(emp_pca_basis,
          labels = paste0("PC", 1:n_components, " (",
                          round(variance_explained * 100), "% of variance)"),
          time = seq(0, 24, by = 0.1),
          title = "Empirical basis from PCA",
          subtitle = "Each component scaled to a maximum absolute value of 1")

