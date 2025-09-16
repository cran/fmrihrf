
#' @export
#' @rdname hrf_from_coefficients
hrf_from_coefficients.HRF <- function(hrf, h, name = NULL, ...) {
  nb <- nbasis(hrf)
  if (length(h) != nb) {
    stop("length(h) must equal nbasis(hrf)")
  }
  weighted_fun <- function(t) {
    vals <- hrf(t)
    if (is.matrix(vals)) {
      drop(vals %*% as.numeric(h))
    } else {
      vals * h[1L]
    }
  }
  if (is.null(name)) {
    name <- paste0(attr(hrf, "name"), "_from_coef")
  }
  as_hrf(
    f      = weighted_fun,
    name   = name,
    nbasis = 1L,
    span   = attr(hrf, "span"),
    params = c(attr(hrf, "params"), list(coefficients = h))
  )
}
