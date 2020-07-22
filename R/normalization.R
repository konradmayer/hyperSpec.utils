#' Normalization
#'
#' @description Perform normalization on objects of class \code{hyperSpec}.
#' @param x A \code{hyperSpec} object
#' @param band numeric or formula; specifying the band or wavelength range to use for band normalization.
#' @param FUN character; "mean" or "sum" to select the function to use for the area normalization.  Using "sum", the integral under the spectrum gets 1, using "mean" the spectrum scales with intensities around 1.
#' @return A \code{hyperSpec} object with normalized intensities
#' @name normalization
NULL

#' @rdname normalization
#' @export
minmax_normalization <- function(x) {
  tmp <- hyperSpec::sweep(x, 1, min, "-")
  out <- hyperSpec::sweep(tmp, 1, max, "/")
  out
}

#' @rdname normalization
#' @export
snv_normalization <- function(x) {
  factors1 <- hyperSpec::apply(x, 1, mean)
  tmp <- hyperSpec::sweep(x, 1, factors1, "-")
  factors2 <- hyperSpec::apply(x, 1, sd)
  out <- hyperSpec::sweep(tmp, 1, factors2, "/")
}

#' @rdname normalization
#' @export
vector_normalization <- function(x) {
  factors1 <- hyperSpec::apply(x, 1, mean)
  tmp <- hyperSpec::sweep(x, 1, factors1, "-")
  factors2 <- hyperSpec::apply(tmp, 1, function(.x) sqrt(sum(.x^2)))
  out <- hyperSpec::sweep(tmp, 1, factors2, "/")
}

#' @rdname normalization
#' @export
area_normalization <- function(x, FUN = c("mean", "sum")) {
  FUN <- match.fun(match.arg(FUN))
  hyperSpec::sweep(x, 1, FUN, "/")
}

#' @rdname normalization
#' @export
band_normalization <- function(x, band) {
  stopifnot(rlang::is_formula(band) | (is.numeric(band) & (length(band) == 1)))
  factors <- hyperSpec::apply(x[, , band], 1, mean)
  hyperSpec::sweep(x, 1, factors, "/")
}
