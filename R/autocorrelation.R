#' Spatial autcorrelation
#'
#' @description This function is a wrapper around \link[ncf]{lisa.nc} to
#'   calculate spatial autocorrelation of a spectral map.
#' @param x a hyperSpec object.
#' @param neigh neighborhood size as defined in \link[ncf]{lisa.nc}.
#' @param colname name of the data column to be added to \code{x@data},
#'   containing the mean correlation within the neighborhood.
#' @param ... additional arguments handed further to \link[ncf]{lisa.nc}.
#'
#' @return a hyperSpec object such as given for \code{x} with an additional
#'   column in \code{x@data}.
#' @export
spatial_ac <- function(x, neigh = 2, colname = "spatial_ac", ...) {
  stopifnot(is_hyperSpecMap(x))
  lisa_out <- ncf::lisa.nc(x$x, x$y, x[[]], neigh = neigh, ...)
  x[, colname] <- lisa_out[["correlation"]]
  x
}

#' Spectral autocorrelation
#'
#' @description Autocorrelation of a spectrum. Generally, a high autocorrelation
#'   is associated with high signal, while white noise is not autocorrelated.
#' @param x a hyperSpec object.
#' @param lag lag at which the autocorrelation is calculated.
#' @param colname name of the data column to be added to \code{x@data},
#'   containing the autocorrelation coefficient.
#'
#' @return a hyperSpec object such as given for \code{x} with an additional
#'   column in \code{x@data}.
#' @export
#'
#' @examples
spectral_ac <- function(x, lag = 1, colname = "spectral_ac") {
  stopifnot(is_hyperSpec(x))
  x[, colname] <- apply(x[[]], 1, function(.x) {
    stats::acf(.x, plot = F)$acf[[lag + 1]]
  })
  x
}
