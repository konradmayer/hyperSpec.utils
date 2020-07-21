#' hyperSpec map to array
#'
#' @description Convert spatially resolved object of class hyperSpec to an array.
#' @param x A \code{hyperSpec} object with variables "x" and "y" in \code{@data}.
#' @return An array with the dimensions "x" "y" and "nwl", with the wavelength vector attached as an attribute.
#' @export
#'
spcmap2array <- function(x) {
  stopifnot(is.hyperSpecMap(x))
  dimensions <- c(x@data[, c("x", "y")] %>% purrr::map(n_distinct), dim(x)["nwl"])
  out <- array(as.vector(x[[]]), dim = dimensions)
  attr(out, "wavelength") <- x@wavelength
  out
}
