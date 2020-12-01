#' Linear combination of spectra
#'
#' @description Use linear regression resp. non negative
#'   linear regression (package {nnls}) without intercept to estimate
#'   coefficients for a linear combination of reference spectra to represent
#'   one or multiple target spectra.
#' @param targets a hyperSpec object.
#' @param references a hyperSpec object containing the spectra to be combined in
#'   in a linear manner.
#' @param references_ids either a single string specifying the column with names
#'   for the references in the object `references`, or a character vector with
#'   a name per reference spectrum.
#' @param nonnegative logical; set TRUE to restrict coefficients to be positive.
#'
#' @return a list containing the following members:
#'   \describe{
#'     \item{data}{the target hyperSpec object with additional columns in `@data` for the coefficients for each of the references as well as summary statistics (RMSE and R2)}
#'     \item{references}{the references object as supplied}
#'     \item{fit}{the return value of the multivariate call to `stats::lm()` or a list of return values of the outputs of `nnls::nnls()` for each target spectrum.}
#'   }
#' @export
linear_combination <- function(targets, references, references_ids = NULL,
                               nonnegative = FALSE) {
  if (is.null(references_ids)) {
    references_ids <- paste0("x", seq_len(nrow(references)))
  } else if (length(references_ids) == 1 & is.numeric(references_ids)) {
    references_ids <- references@data[, references_ids]
  } else if (length(references_ids) == 1 & is.character(references_ids)) {
    references_ids <- references_ids
  } else {
    stopifnot(length(references_ids) == nrow(references))
  }


  y <- t(targets[[]])
  x <- as.data.frame(t(references[[]]))
  colnames(x) <- references_ids

  if (nonnegative) {
    # lm has no option to restrict coefs to positive values. As coefs can be
    # interpreted as quantity of a chemical compound in a mixture, only positive
    # values are meaningful. The argument nonnegative allows for this restriction.
    # init
    coefs <- as.data.frame(matrix(ncol = nrow(references), nrow = nrow(targets)))
    colnames(coefs) <- paste("coef", references_ids, sep = "_")
    rmse <- r2 <- numeric(nrow(targets))
    fit <- list()
    # no multivariate fit (as for lm below) possible, therefore looped
    for (i in seq_len(ncol(y))) {
      fit_tmp <- nnls::nnls(as.matrix(x), y[, i])
      coefs[i, ] <- fit_tmp$x # coefs under member "x"
      rmse[[i]] <- sqrt(colMeans(fit_tmp$residuals^2))
      r2[[i]] <- 1 - (sum(fit_tmp$residuals^2) / sum((y[, i] - mean(y[, i]))^2))
      fit[[i]] <- fit_tmp
    }
  } else {
    fit <- stats::lm(y ~ 0 + ., data = x)
    coefs <- t(stats::coef(fit))
    colnames(coefs) <- paste("coef", colnames(coefs), sep = "_")
    rmse <- sqrt(colMeans(fit$residuals^2))
    r2 <- vapply(summary(fit), function(x) x[["r.squared"]], numeric(1))
  }

  summary_stats <- data.frame(rmse = rmse, r2 = r2)
  additional_cols <- cbind(coefs, summary_stats)
  colnames(additional_cols) <- paste("LC", colnames(additional_cols), sep = "_")


  targets@data <- cbind(targets@data, additional_cols)

  out <- list(
    data = targets,
    references = references,
    fit = fit
  )

  class(out) <- "lc"
  out
}
