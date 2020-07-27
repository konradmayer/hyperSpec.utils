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
