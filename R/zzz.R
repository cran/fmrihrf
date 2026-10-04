.onLoad <- function(libname, pkgname) {
  # knitr is suggested, not imported: register the knit_print method for
  # fmrihrf_plot whenever knitr is (or later becomes) loaded.
  register <- function(...) {
    registerS3method("knit_print", "fmrihrf_plot", knit_print.fmrihrf_plot,
                     envir = asNamespace("knitr"))
  }
  if (isNamespaceLoaded("knitr")) {
    register()
  } else {
    setHook(packageEvent("knitr", "onLoad"), register)
  }
  invisible()
}
