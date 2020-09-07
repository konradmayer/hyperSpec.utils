#' Dimensions of a hyperSpec map
#'
#' @param x hyperSpec object with variables "x" and "y" present in \code{@data}.
#'
#' @return a named vector with the number of x and y positions as well as the number of descrete wavelengths ("nwl").
#' @export
#'
#' @examples
#' library(hyperSpec)
#' spcmap_dim(chondro)
spcmap_dim <- function(x) {
  if (!is_hyperSpecMap(x)) {
    stop("This function is only applicable to spatially resolved hyperSpec objects")
  }
  c(
    vapply(x@data[, c("x", "y")], function(x) length(unique(x)), numeric(1)),
    dim(x)["nwl"]
  )
}


#' Aggregate wavelengths
#'
#' @description select certain bands of aggregate wavelength ranges of a hyperSpec object.
#' @param x a hyperSpec object.
#' @param FUN function to be used for aggregation.
#' @param ... (named) selection of bands or band ranges.
#'
#' @return a hyperspec object with aggregated spc
#' @export
#'
#' @examples
#' \dontrun{
#' aggregate_wl(chondro, "mean", band_a = 938, band_b = 1230 ~ 1260)
#' }
aggregate_wl <- function(x, FUN = "mean", ...) {
  aggregates <- list(...)
  FUN <- match.fun(FUN)
  # replace lapply with vapply here?
  tmp <- lapply(aggregates, function(.x) hyperSpec::apply(x[, , .x], 1, FUN))
  out <- do.call(cbind, tmp)
  colnames(out$spc) <- names(aggregates)
  out
}



#' Plotting Spectra
#'
#' @description Simple wrapper for \code{\link[hyperSpec]{plotspc}} to execute
#'   the function with reversed x-axis as default.
#' @inheritParams hyperSpec::plotspc
#' @param ... additional arguments handed to \code{\link[hyperSpec]{plotspc}}.
#' @inherit hyperSpec::plotspc return
#' @author NULL
#' @export

plotspc_rev <- function(object, ...) {
  hyperSpec::plotspc(object, wl.reverse = TRUE, ...)
}



#' Test if the wavelength vector of a hyperSpec object is equally spaced
#'
#' @param x a hyperSpec object
#'
#' @return logical
#' @export
is_wl_equidistant <- function(x) {
  all(diff(sort(hyperSpec::wl(x)), differences = 2) == 0)
}



#' Test if the wavelength vector of a hyperSpec object is in order (ascending)
#'
#' @param x a hyperSpec object
#'
#' @return logical
#' @export
is_wl_ordered <- function(x) {
  wl <- hyperSpec::wl(x)
  all(wl == sort(wl))
}
