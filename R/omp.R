
#' Spectral mixture analysis: Orthogonal Matching Pursuit (OMP)
#'
#' @description Use the Orthogonal Matching Pursuit to estimate sparse
#'   coefficients for a linear combination of reference spectra to represent
#'   one or multiple target spectra.
#' @param targets a hyperSpec object.
#' @param references a hyperSpec object containing the spectra to be combined in
#'   in a linear manner.
#' @param references_ids either a single string or numeric specifying the column with names
#'   for the references in the object `references`, or a character vector with
#'   a name per reference spectrum.
#' @param tol The tolerance value to terminate the algorithm, see
#'   \code{\link[Rfast]{ompr}}. As the function is run with
#'   \code{method = "SSE"}, the algorithm stops if the relative change in the
#'   sum of squared errors is equal to or less than \code{tol}.
#' @param parallel logical; If this is \code{TRUE},
#'   \code{\link[parallel]{parApply}} is used. Please be aware that the function
#'   handles starting and stoping of the cluster (type = "PSOCK").
#' @param ncores number of logical cores to use for parallel processing.
#'   Be aware that more is not necessarily better, as starting new processes
#'   adds overhead which potentially is bigger than the time savings from
#'   parallel processing.
#' @details This is a wrapper around \code{\link[Rfast]{ompr}}. The function is
#'   run with \code{method = "SSE"}.
#' @references Pati, Y. C., R. Rezaiifar, Y. C. Pati R. Rezaiifar, and P. S. Krishnaprasad. “Orthogonal Matching Pursuit: Recursive Function Approximation with Applications to Wavelet Decomposition.” In Proceedings of the 27 Th Annual Asilomar Conference on Signals, Systems, and Computers, 40–44, 1993.
#' @return a list containing the following members:
#'   \describe{
#'     \item{coefficients}{coefficient matrix}
#'     \item{basis}{the references object as supplied}
#'     \item{fit}{the return value of the multivariate call to `stats::lm()` or a list of return values of the outputs of `nnls::nnls()` for each target spectrum.}
#'   }
#' @export

omp <- function(targets, references, tol = 0.2, references_ids = NULL,
                parallel = FALSE, ncores = 4) {
  references_ids <- get_references_ids(references, references_ids)

  if (parallel) {
    refmat <- t(references[[]])
    parapply_omp_fun <- function(x, refmat, tol) {
      ompr_out <- Rfast::ompr(x, refmat,
        method = "sse", xstand = TRUE,
        ystand = TRUE, tol = tol
      ) # while the RSS dont get better than specified fraction (tol)


      out <- numeric(ncol(refmat))
      if (nrow(ompr_out$info) > 1) {
        omp_coefs <- nnls::nnls(refmat[, ompr_out$info[-1, 1], drop = FALSE], x)$x # coefs in member x
        out[ompr_out$info[-1, 1]] <- as.numeric(omp_coefs)
      }

      out
    }

    clust <- parallel::makeCluster(ncores)
    print(clust)
    on.exit(parallel::stopCluster(clust))
    coefs_outmat <- t(parallel::parApply(clust, targets[[]], 1, parapply_omp_fun, refmat, tol))
  } else {
    ntargets <- nrow(targets)

    coefs_outmat <- matrix(0, ncol = nrow(references), nrow = ntargets) # init out object
    x <- t(references[[]])
    for (i in seq_len(ntargets)) {
      y <- as.numeric(t(targets[i, ][[]]))
      ompr_out <- Rfast::ompr(y, x,
        method = "sse", ystand = TRUE,
        xstand = TRUE, tol = tol
      ) # while the RSS dont get better than specified fraction (tol)

      if (nrow(ompr_out$info) > 1) {
        refidx <- ompr_out$info[-1, 1]
        omp_coefs <- nnls::nnls(x[, refidx, drop = FALSE], y)[["x"]] # coefs in member x
        coefs_outmat[i, refidx] <- as.numeric(omp_coefs)
      }
    }
  }

  colnames(coefs_outmat) <- paste("coef", references_ids, sep = "_")

  out <- list(
    coefficients = coefs_outmat,
    basis = references
  )

  class(out) <- "omp"
  out
}
