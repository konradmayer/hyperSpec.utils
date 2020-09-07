
# S3 methods --------------------------------------------------------------

# obsolete - see S4 method for diff below
# #' differences on hyperSpec spectra (S3 method)
# #'
# #' @description S3 method for the generic function \code{diff()} and the hyperSpec class.
# #'
# #' @param x a hyperSpec object.
# #' @details All columns of \code{x@data$spc} are kept and spectra padded with NAs accordingly.
# #' @inherit base::diff
# #' @return a hyperSpec object.
# #' @export
# #'
# #' @examples \dontrun{
# #' plotspc(diff(chondro, 1, 1))
# #' }
# #'
# diff.hyperSpec <- function(x, lag = 1, differences = 1) {
#   hyperSpec::apply(x, 1, function(.x) c(rep(NA, times = differences * lag),
#                                         diff(.x, lag = lag, differences = differences)))
# }



# S4 methods --------------------------------------------------------------

#' differences on hyperSpec spectra (S4 method)
#'
#' @description S4 method for the generic function \code{diff()} and the hyperSpec class.
#' @details All columns of \code{x@data$spc} are kept and spectra padded with NAs accordingly.
#' @inheritParams base::diff
#' @return a hyperSpec object.
#' @aliases diff
#' @export
#' @seealso \code{\link[base]{diff}}
#' @examples
#' \dontrun{
#' plotspc(diff(chondro, 1, 1))
#' }
setMethod(
  f = "diff",
  signature = "hyperSpec",
  definition = function(x, lag = 1, differences = 1) {
    x@data[["spc"]] <- t(apply(x@data$spc, 1, function(.x) {
      c(
        rep(NA, times = differences * lag),
        diff(.x, lag = lag, differences = differences)
      )
    }))
    hyperSpec::chk.hy(x)
    return(x)
  }
)
