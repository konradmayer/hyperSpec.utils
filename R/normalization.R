minmax_normalization <- function(x) {
  tmp <- hyperSpec::sweep(x, 1, min, "-")
  out <- hyperSpec::sweep(x, 1, max, "/")
  out
}

snv_normalization <- function(x) {
  factors1 <- hyperSpec::apply(x, 1, mean)
  tmp <- hyperSpec::sweep(x, 1, factors1, "-")
  factors2 <- hyperSpec::apply(x, 1, sd)
  out <- hyperSpec::sweep(tmp, 1, factors2, "/")
}


vector_normalization <- function(x) {
  factors1 <- hyperSpec::apply(x, 1, mean)
  tmp <- hyperSpec::sweep(x, 1, factors1, "-")
  factors2 <- hyperSpec::apply(tmp, 1, function(.x) sqrt(sum(.x^2)))
  out <- hyperSpec::sweep(tmp, 1, factors2, "/")
}


area_normalization <- function(x, FUN = c("mean", "sum")) {
  FUN <- match.fun(match.arg(FUN))
  hyperSpec::sweep(x, 1, FUN, "/")
}


band_normalization <- function(x, band) {
  stopifnot(is_formula(band) | (is.numeric(band) & (length(band) == 1)))
  factors <- hyperSpec::apply(x[, , band], 1, mean)
  hyperSpec::sweep(x, 1, factors, "/")
}
