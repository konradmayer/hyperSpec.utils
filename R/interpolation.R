
# helpers are unexported functions from hyperSpec

`.wl<-` <- function(x, value) {
  x@wavelength <- value
  spc <- .fix_spc_colnames(x)
  x
}

.fix_spc_colnames <- function(spc) {
  colnames(spc@data$spc) <- signif(spc@wavelength, digits = 6)
  spc
}


#' Interpolation of spectra onto a new wavelength axis
#'
#' @param spc an object of class \code{hyperSpec}
#' @param newx numeric; a new wavelength (or wavenumber) vector
#' @param ... additional arguments to be handed further to \link[hyperSpec]{spc.smooth.spline} or \link[stats]{approx}.
#' @name interpolation
NULL

#' @rdname interpolation
#' @export
spc_linear_interpolation <- function(spc, newx = hyperSpec::wl(spc), ...) {
  .linear_interpolation <- function(x, y, newx) {
    pts <- !is.na(y)
    stats::approx(x[pts], y[pts], newx, ...)$y
  }
  spc <- hyperSpec::orderwl(spc)
  newspc <- matrix(NA_real_, ncol = length(newx), nrow = nrow(spc))
  i <- rowSums(is.na(spc@data$spc)) < hyperSpec::nwl(spc)
  newspc[i, ] <- t(apply(spc@data$spc[i, , drop = FALSE],
    1, .linear_interpolation,
    x = spc@wavelength, newx = newx
  ))
  if (any(is.na(newspc[i, ]))) {
    warning("NAs generated. Probably newx was outside the spectral range covered by spc.")
  }
  spc@data$spc <- newspc
  .wl(spc) <- newx
  methods::validObject(spc)
  spc
}

#' @rdname interpolation
#' @export
spc_spline_interpolation <- function(spc, newx = hyperSpec::wl(spc), ...) {
  hyperSpec::spc.smooth.spline(spc, newx, spar = 0, all.knots = TRUE, ...)
}
