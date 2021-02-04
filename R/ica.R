
#' Spectral unmixing: Independent Component Analysis (ICA)
#'
#' @description This is a wrapper around \code{\link[fastICA]{fastICA}}
#'   from the package {fastICA} to be used with objects of class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer; number of pure components.
#' @param prefix character; a prefix to name the pure spectra.
#' @param ... further arguments to be passed down to \code{\link[fastICA]{fastICA}}.
#' @references A. Hyvarinen, "Fast and robust fixed-point algorithms for independent component analysis," in IEEE Transactions on Neural Networks, vol. 10, no. 3, pp. 626-634, May 1999 \url{https://doi.org/10.1109/72.761722}
#' @return a list with the following components:
#'   \describe{
#'     \item{coefficients}{coefficient matrix}
#'     \item{basis}{a hyperSpec object containing the basis (component) spectra}
#'     \item{fit}{the ica fit as returned by  \code{\link[fastICA]{fastICA}}}
#'   }
#' @seealso \code{\link[fastICA]{fastICA}}
#' @export

ica <- function(x, ncomp, prefix = "basis", ...) {
  results <- fastICA::fastICA(x[[]], ncomp)

  component_nm <- paste0("NMF_", prefix, seq_len(ncomp))
  components <- methods::new("hyperSpec",
    spc = results$A, wavelength = hyperSpec::wl(x),
    data = data.frame(basis = component_nm)
  )
  params <- results$S
  colnames(params) <- component_nm

  out <- list(
    coefficients = params,
    basis = components,
    fit = results
  )
  class(out) <- "ica"
  out
}
