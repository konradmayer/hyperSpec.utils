
# helpers -----------------------------------------------------------------

# from package alsace, but all arguments made available in function call and naming changed (prefix)
.doALS_custom <- function(Xl, PureS, maxiter = 100,
                          normS = 0.5,
                          nonnegS = TRUE,
                          optS1st = FALSE,
                          nonnegC = TRUE,
                          uniC = FALSE,
                          baseline = FALSE,
                          prefix = "component") {
  Cini <- lapply(Xl, function(xl) xl[, 1:ncol(PureS)])
  utils::capture.output(result <- ALS::als(
    PsiList = Xl,
    CList = Cini,
    S = PureS,
    maxiter = maxiter,
    normS = normS,
    nonnegS = nonnegS,
    optS1st = optS1st,
    nonnegC = nonnegC,
    uniC = uniC,
    baseline = baseline
  ))
  colnames(result$S) <- paste0(prefix, 1:ncol(PureS))
  for (i in 1:length(result$CList)) colnames(result$CList[[i]]) <- colnames(result$S)
  predicted.values <- lapply(1:length(Xl), function(ii) {
    Xl[[ii]] -
      result$resid[[ii]]
  })
  predvals2 <- sum(unlist(predicted.values)^2)
  npoints <- prod(c(
    length(result$CList), nrow(result$CList[[1]]),
    nrow(result$S)
  ))
  result$summ.stats <- list(
    lof = 100 * sqrt(result$rss / predvals2),
    rms = sqrt(result$rss / npoints), r2 = 1 - result$rss / predvals2
  )
  class(result) <- "ALS"
  result
}


# exported functions ------------------------------------------------------

#' Orthogonal Projection Approach
#'
#' @description find the most dissimilar spectra in a hyperSpec object. This
#'   function is a customized version of the function \link[alsace]{opa} from
#'   the package {alsace} to work with objects of the class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer, number of most dissimilar spectra to return.
#' @param return_scaled if TRUE, unit-length scaled spectra are returned,
#'   if FALSE the original intensities are retained.
#'
#' @return a hyperSpec object containing ncomp spectra.
#' @export
#'
opa <- function(x, ncomp = NULL, return_scaled = TRUE) {
  x_mat <- x[[]]
  lambdas <- colnames(x_mat)
  x_mat <- t(apply(x_mat, 1, function(xx) {
    xx / rep(
      sqrt(crossprod(xx)),
      length(xx)
    )
  }))
  Xref <- matrix(0, ncomp, ncol(x_mat))
  huhn <- colMeans(x_mat)
  Xref[1, ] <- huhn / rep(sqrt(crossprod(huhn)), length(huhn))
  ids <- integer(length = ncomp)
  for (i in seq_len(ncomp)) {
    Xs <- lapply(1:nrow(x_mat), function(ii, xx, xref) {
      rbind(
        xref,
        xx[ii, ]
      )
    }, x_mat, Xref[1:(i - 1), ])
    dissims <- sapply(Xs, function(xx) det(tcrossprod(xx)))
    id <- which.max(dissims)
    ids[[i]] <- id
    Xref[i, ] <- x_mat[id, ]
  }

  if (return_scaled) {
    # return unit-length scaled - such as function alsace::opa scales spectra
    hyperSpec::apply(x[ids, ], 1, function(xx) xx / rep(sqrt(crossprod(xx)), length(xx)))
  } else {
    # return original values
    x[ids, ]
  }
}



#' Spectral unmixing: Alternating least squares multivariate curve resolution (MCR-ALS)
#'
#' @description This is a wrapper around a modified version of \link[alsace]{doALS}
#'   from the package {alsace} to be used with objects of class hyperSpec.
#' @param x a hyperSpec object.
#' @param ncomp integer; number of pure components, does not need to be specified
#'   if S_init is provided.
#' @param S_init a hyperSpec object; initial estimates of pure spectral components,
#'   if S_init is NULL, initial components are estimated
#'   using \code{\link{opa}}.
#' @param prefix character; a prefix to name the pure spectra.
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{data}{the input hyperSpec object with additional "component" columns in @data, holding the coefficients}
#'     \item{component_spectra}{a hyperSpec object containing the component spectra}
#'     \item{summary_stats}{a named vector containing lof, rms and r2}
#'     \item{fit}{the als fit as returned by the custom wrapper around `ALS::als()`}
#'   }
#' @seealso \link[alsace]{doALS}, \link[ALS]{als}
#' @export
#'

als <- function(x, ncomp = NULL, S_init = NULL, prefix = "component") {
  if (!is_hyperSpec(x)) {
    stop("x needs to be an object of class hyperSpec.")
  }

  if (!is.null(S_init)) {
    if ((!is.null(ncomp)) && (ncol(S_init) != ncomp)) {
      stop("provided ncomp does not match the number of colums in the provided S_init")
    }
  } else { # if S_init not provided, start values are calculated with opa (changed from alsace::opa)
    S_init <- opa(x, ncomp = ncomp)
  }

  fit <- .doALS_custom(list(x[[]]), t(S_init[[]]), prefix = prefix)

  params <- fit$CList[[1]]
  colnames(params) <- paste("ALS", colnames(params), sep = "_")
  x@data <- cbind(x@data, params)
  components <- methods::new("hyperSpec",
    spc = t(fit$S),
    wavelength = hyperSpec::wl(x),
    data = data.frame(
      component =
        paste0("ALS_", prefix, seq_len(ncomp))
    )
  )

  out <- list(
    data = x,
    component_spectra = components,
    summary_stats = unlist(fit$summ.stats),
    fit = fit
  )
  class(out) <- "als"
  out
}

summary.als <- function(x) {
  cat(
    "ALS fit with", nrow(x$data), "samples with",
    nrow(x$component_spectra), "components.", "\nEach hyperspec object contains",
    hyperSpec::nwl(x$data), "wavelengths.\n"
  )
  cat(
    "\tRMSE of fit:", round(x$summary_stats[["rms"]], 5),
    "\n"
  )
  cat("\tLOF: ", round(x$summary_stats[["lof"]], 2), "%\n", sep = "")
  cat("\tR2: ", round(x$summary_stats[["r2"]], 5), "\n", sep = "")
  invisible()
}






# out <- als(dat, 4)
# summary(out)
# plotspc(out$component_spectra, stacked = TRUE)
# plotmap_viridis(out$data, component4 ~ x * y)
