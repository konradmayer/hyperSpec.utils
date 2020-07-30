#' cubeview for hyperSpec objects
#'
#' @description This function is a wrapper for the function
#'   \href{https://github.com/r-spatial/cubeview}{cubeview} to be used with
#'   objects of the class \code{hyperSpec}. This function allows you to interactively
#'   explore a spectral hypercube. Use the LEFT / RIGHT arrow keys and DOWN / UP
#'   arrow keys to set the slicing positions of the x-axis and y-axis,
#'   PAGE_DOWN / PAGE_UP keys for the z-axis. Press SPACE to show or hide guides.
#'   Use the left mouse-button to rotate, right mouse button to pan the hypercube.
#'   If you use the function inside RStudio, the Viewer pane may not display the cube.
#'   In that case, click on "Show in new window" to open it in a browser.
#' @param x a hyperSpec object.
#' @param ... additional arguments passed on to \code{cubeview::cubeview()}.
#' @export
#'
#' @examples
#' \dontrun{
#' library(hyperSpec)
#' library(hyperSpec.utils)
#' cubeview(chondro)
#' }
cubeview <- function(x, ...) {
  cubeview::cubeview(raster::brick(spcmap2array(x)), ...)
}
