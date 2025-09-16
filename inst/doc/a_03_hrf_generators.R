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
library(ggplot2)
library(dplyr)
library(tidyr)

## ----list-generators----------------------------------------------------------
list_available_hrfs(details = TRUE) %>%
  dplyr::filter(type == "generator")

## ----create-basis-------------------------------------------------------------
# Create a B-spline basis using gen_hrf
bs8 <- gen_hrf(hrf_bspline, N = 8, span = 32)
print(bs8)

## ----eval-basis---------------------------------------------------------------
times <- seq(0, 32, by = 0.5)
mat <- bs8(times)
head(mat)

## ----fir-basis----------------------------------------------------------------
# Use the pre-defined FIR basis or create one with gen_hrf
fir10 <- HRF_FIR  # Pre-defined FIR with 12 basis functions
resp <- fir10(times)

fir_df <- data.frame(Time = times, resp)
fir_long <- tidyr::pivot_longer(fir_df, -Time)

ggplot(fir_long, aes(Time, value, colour = name)) +
  geom_line(linewidth = 1) +
  labs(title = "Finite Impulse Response Basis",
       x = "Time (s)", y = "Response") +
  theme_minimal() +
  theme(legend.position = "none")

## ----gethrf, eval=FALSE-------------------------------------------------------
# # Internal usage only:
# # custom_fir <- getHRF("fir", nbasis = 6, span = 18)
# # custom_fir
# ``` -->
# 
# ## Summary
# 
# Generator functions are simple factories that let you customise flexible HRF
# bases. They return normal `HRF` objects, which means you can evaluate them,
# combine them with decorators, or insert them into regressors just like the
# built-in HRFs.

