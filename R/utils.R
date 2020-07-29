#' Dimensions of a hyperSpec map
#'
#' @param x hyperSpec object with variables "x" and "y" present in \code{@data}.
#'
#' @return a named vector with the number of x and y positions as well as the number of descrete wavelengths ("nwl")
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
  out <- hyperSpec:::cbind.hyperSpec(tmp)
  colnames(out$spc) <- names(aggregates)
  out
}
