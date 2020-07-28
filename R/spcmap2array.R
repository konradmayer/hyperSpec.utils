#' hyperSpec map to array
#'
#' @description Convert spatially resolved object of class hyperSpec to an array.
#' @param x A \code{hyperSpec} object with variables "x" and "y" in \code{@data}.
#' @param set_dimnames logical; specify whether values for x, y and wavelength should be used as dimnames.
#' @param sort logical; specify whether dimensions x (increasing), y (decreasing) and wavelengths (increasing) should be sorted prior conversion.
#' @return An array with the dimensions "x" "y" and "nwl", with the wavelength vector attached as an attribute.
#' @export
#'
spcmap2array <- function(x, set_dimnames = FALSE, sort = FALSE) {
  # input validation
  stopifnot(is.hyperSpecMap(x))

  # sort
  if (sort) {
    x <- x[order(x$x, decreasing = FALSE), ]
    x <- x[order(x$y, decreasing = TRUE), ]
    x <- hyperSpec::orderwl(x)
  }
  # create array
  dimensions <- spcmap_dim(x)
  out <- array(as.vector(x[[]]), dim = dimensions)
  attr(out, "wavelength") <- x@wavelength
  if (set_dimnames) {
    dimnames(out) <- list(
      x = unique(x$x), y = unique(x$y),
      nwl = x@wavelength
    )
  }

  out
}
