#' Spectral unmixing: Nonnegative matrix factorization (NMF)
#'
#' @description This is a wrapper around \link[NMF]{nmf}
#'   from the package {NMF} to be used with objects of class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer; number of pure components (factorization rank).
#' @param nrun number of runs
#' @param method character; use `NMF::getNMFMethod()` to get valid Methods and
#'   refer to the documentation of `NMF::nmf()` for further information.
#' @param seed character; "opa" (Orthogonal Projection Approach, see \code{\link{omp}}) or one of `NMF::getNMFSeed()`. Refer to the documentation
#'   of `NMF::nmf()` for further information.
#' @param prefix character; a prefix to name the pure spectra.
#' @param ... further arguments to be passed down to \link[NMF]{nmf}.
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{coefficients}{coefficient matrix}
#'     \item{basis}{a hyperSpec object containing the basis (component) spectra}
#'     \item{fit}{the nmf fit as returned by `NMF::nmf()`}
#'   }
#' @seealso \link[NMF]{nmf}
#' @export

nmf <- function(x, ncomp = NULL, nrun = 1, method = "brunet", seed = "opa", prefix = "basis", ...) {
  if (!is_hyperSpec(x)) {
    stop("x needs to be an object of class hyperSpec.")
  }

  if (is.null(ncomp)) {
    stop("please specify the number of end members using the argument ncomp")
  }

  # additional seeding method using opa
  if (seed == "opa") {
    opa_basis <- opa(x, ncomp)
    suppressWarnings(lc_opa <- linear_combination(x, opa_basis))
    suppressWarnings(seed <- NMF::nmfModel(rank = ncomp, H = t(lc_opa$coefficients), W = t(lc_opa$basis[[]])))
  }


  fit <- NMF::nmf(t(x[[]]), ncomp, nrun = nrun, method = method, seed = seed, ...)

  component_nm <- paste0("NMF_", prefix, seq_len(ncomp))
  components <- methods::new("hyperSpec",
    spc = t(NMF::basis(fit)), wavelength = hyperSpec::wl(x),
    data = data.frame(basis = component_nm)
  )
  params <- t(NMF::coef(fit))
  colnames(params) <- component_nm


  out <- list(
    coefficients = params,
    basis = components,
    fit = fit
  )
  class(out) <- "nmf"
  out
}
