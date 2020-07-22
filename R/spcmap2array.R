#' hyperSpec map to array
#'
#' @description Convert spatially resolved object of class hyperSpec to an array.
#' @param x A \code{hyperSpec} object with variables "x" and "y" in \code{@data}.
#' @return An array with the dimensions "x" "y" and "nwl", with the wavelength vector attached as an attribute.
#' @export
#'
spcmap2array <- function(x) {
  stopifnot(is.hyperSpecMap(x))
  dimensions <- c(
    vapply(x@data[, c("x", "y")], function(x) length(unique(x)), numeric(1)),
    dim(x)["nwl"]
  )
  out <- array(as.vector(x[[]]), dim = dimensions)
  attr(out, "wavelength") <- x@wavelength
  out
}
