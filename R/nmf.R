#' Spectral unmixing: Nonnegative matrix factorization (NMF)
#'
#' @description This is a wrapper around \link[NMF]{nmf}
#'   from the package {NMF} to be used with objects of class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer; number of pure components (factorization rank).
#' @param nrun number of runs
#' @param method character; use `NMF::getNMFMethod()` to get valid Methods and
#'   refer to the documentation of `NMF::nmf()` for further information.
#' @param seed character; one of `NMF::getNMFSeed()`. Refer to the documentation
#'   of `NMF::nmf()` for further information.
#' @param prefix character; a prefix to name the pure spectra.
#' @param ... further arguments to be passed down to \link[NMF]{nmf}.
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{data}{the input hyperSpec object with additional "component" columns in @data, holding the coefficients}
#'     \item{component_spectra}{a hyperSpec object containing the component spectra}
#'     \item{fit}{the nmf fit as returned by `NMF::nmf()`}
#'   }
#' @seealso \link[NMF]{nmf}
#' @export

nmf <- function(x, ncomp = NULL, nrun = 1, method = "brunet", seed = "ica", prefix = "component", ...) {
  if (!is_hyperSpec(x)) {
    stop("x needs to be an object of class hyperSpec.")
  }

  if (is.null(ncomp)) {
    stop("please specify the number of end members using the argument ncomp")
  }

  fit <- NMF::nmf(x[[]], ncomp, nrun = nrun, method = method, seed = seed, ...)

  # coef is matrix H
  components <- methods::new("hyperSpec",
    spc = NMF::coef(fit), wavelength = hyperSpec::wl(x),
    data = data.frame(component = paste0("NMF_", prefix, seq_len(ncomp)))
  )
  # basis is matrix W
  params <- stats::setNames(as.data.frame(NMF::basis(fit)), paste0("NMF_", prefix, seq_len(ncomp)))
  x@data <- cbind(x@data, params)


  out <- list(
    data = x,
    component_spectra = components,
    fit = fit
  )
  class(out) <- "nmf"
  out
}
