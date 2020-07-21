#' Test if the object is of class hyperSpec
#'
#' @description This function returns TRUE for objects of class \code{hyperSpec}, and FALSE for all other objects.
#' The variant \code{is.hyperSpecMap} returns only TRUE if additionally variables "x" and "y" are present in \code{@data}.
#' @aliases is.hyperSpec is.hyperSpecMap
#' @param x An object
#' @name is_hyperSpec
NULL

#' @rdname is_hyperSpec
#' @export
is_hyperSpec <- is.hyperSpec <- function(x) {
  inherits(x, "hyperSpec")
}

#' @rdname is_hyperSpec
#' @export
is_hyperSpecMap <- is.hyperSpecMap <- function(x) {
  inherits(x, "hyperSpec") & all(c("x", "y") %in% names(x@data))
}
