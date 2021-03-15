#' Spectral unmixing: Vertex component analysis (VCA)
#'
#' @description This is a wrapper around \code{vca} in package \code{{unmixR}}
#'   to be used with objects of class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer; number of pure components.
#' @param prefix character; a prefix to name the pure spectra.
#' @param ... further arguments to be passed down to \code{vca} in package \code{{unmixR}}
#' @details Use \code{method = "05"} for the algorithm described in Nascimento et al.
#' 2005 or \code{method = "Lopez2012"} for the algorithm described in Lopez et al.
#' 2012. Refer to the documentation of \code{vca} in package \code{{unmixR}}
#' for further details.
#'
#' @references Nascimento, J.M.P. and Bioucas Dias, J.M. "Vertex component
#'   analysis: a fast algorithm to unmix hyperspectral data," Geoscience and
#'   Remote Sensing, vol. 43, no. 4, pp. 898-910, April 2005,
#'   doi: 10.1109/TGRS.2005.844293
#'
#'   Lopez, S., Horstrand, P., Callico, G.M., Lopez J.F. and
#'   Sarmiento, R., "A Low-Computational-Complexity Algorithm for
#'   Hyperspectral Endmember Extraction: Modified Vertex Component Analysis,"
#'   Geoscience & Remote Sensing Letters, IEEE, vol. 9 no. 3 pp. 502-506, May 2012
#'   doi: 10.1109/LGRS.2011.2172771
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{coefficients}{coefficient matrix}
#'     \item{basis}{a hyperSpec object containing the basis (component) spectra}
#'   }
#' @export
#'

vca <- function(x, ncomp, prefix = "basis", ...) {
  out_vcs <- unmixR::vca(x, ncomp, ...)
  component_nm <- paste0("VCA_", prefix, seq_len(ncomp))
  basis <- methods::new("hyperSpec",
    spc = x[out_vcs$indices, ][[]],
    wavelength = hyperSpec::wl(x),
    data = data.frame(basis = component_nm)
  )

  out <- linear_combination(x, basis)[1:2]
  colnames(out$coefficients) <- component_nm
  class(out) <- "vca"
  out
}
