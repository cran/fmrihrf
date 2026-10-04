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
library(ggplot2)
library(dplyr)
library(tidyr)

## ----list-generators---------------------------------------------
list_available_hrfs(details = TRUE) %>%
  dplyr::filter(type == "generator")

## ----create-basis------------------------------------------------
# Create a B-spline basis using gen_hrf
bs8 <- gen_hrf(hrf_bspline, N = 8, span = 32)
print(bs8)

## ----eval-basis--------------------------------------------------
times <- seq(0, 32, by = 0.5)
mat <- bs8(times)
head(mat)

## ----fir-basis, fig.alt="Ten FIR basis functions, each a 2 second boxcar of height 1 labelled B1 to B10; bin k covers 2(k-1) to 2k seconds after the event, so together they tile 0 to 20 seconds."----
fir10 <- hrf_fir_generator(nbasis = 10, span = 20)
print(fir10)

plot_hrfs(fir10, time = seq(0, 22, by = 0.02),
          title = "FIR basis: 10 bins of 2 s")

## ----gethrf, eval=FALSE------------------------------------------
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

