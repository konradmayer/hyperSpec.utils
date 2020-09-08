
# helpers -----------------------------------------------------------------


#' Trapezoid integration
#'
#' @param x numeric.
#' @param y numeric.
#'
#' @return numeric; the approximated integral.
#' @export
integrate_trapezoid <- function(x, y) {
  # from pracma::trapz
  if (missing(y)) {
    if (length(x) == 0) {
      return(0)
    }
    y <- x
    x <- seq(along = x)
  }
  if (length(x) == 0 && length(y) == 0) {
    return(0)
  }
  if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) ||
    is.complex(y))) {
    stop("Arguments 'x' and 'y' must be real or complex vectors.")
  }
  m <- length(x)
  if (length(y) != m) {
    stop("Arguments 'x', 'y' must be vectors of the same length.")
  }
  if (m <= 1) {
    return(0)
  }
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2 * m
  p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
  p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
  return(0.5 * (p1 - p2))
}

#' Rectangular (midpoint) integration
#'
#' @param x numeric.
#' @param y numeric.
#' @param correct_edges logical; if set TRUE, the widths of the segments at the
#'   range edges (the first and last point) are corrected to only incorporate
#'   the half width to the right, resp. to the left of the point.
#'
#' @return numeric; the approximated integral.
#' @export
integrate_rectangular <- function(x, y, correct_edges = TRUE) # if correct edges = TRUE, result same as trapezoidal integration
{
  if (missing(y)) {
    if (length(x) == 0) {
      return(0)
    }
    y <- x
    x <- seq(along = x)
  }
  if (length(x) == 0 && length(y) == 0) {
    return(0)
  }
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Arguments 'x' and 'y' must be real vectors.")
  }
  m <- length(x)
  if (length(y) != m) {
    stop("Arguments 'x', 'y' must be vectors of the same length.")
  }
  if (m <= 1) {
    return(0)
  }

  # define x width belonging to a y value
  halfdiffs <- diff(x) / 2
  left <- c(halfdiffs[1], halfdiffs)
  right <- c(halfdiffs, halfdiffs[length(halfdiffs)])
  widths <- left + right
  # correct the leftmost and rightmost value to not extend the x range by halfwidth
  if (correct_edges) {
    widths[c(1, length(widths))] <- halfdiffs[c(1, length(halfdiffs))]
  }
  # calculate integral
  integral <- sum(widths * y)
  integral
}


# alternatives derived from https://www.rdocumentation.org/packages/fda.usc/versions/2.0.2/source
# integrate_trapezoid <- function(x, y) {
#   n <- length(x)
#   idx <- 2:n
#   value<-as.double((x[idx]-x[idx-1])%*%(y[idx]+y[idx-1]))/2
#   value
# }
# integrate_simpson_composite <- function(x, y, n = length(x)) {
#   n <- 2*n-1
#   ## use linear approximation to get equally spaced x values
#   app <- approx(x,y,n=n)
#   x <- app$x
#   y <- app$y
#   h <- (max(x)-min(x))/(n-1)
#   value <- (h/3) * (y[n] + y[1] + 2*sum(y[2*(1:((n-1)/2))+1]) + 4*sum(y[2*(1:((n-1)/2))]))
#   value
# }
# integrate_simpson_extended <- function(x, y, n = length(x)) {
#   n <- 2*n-1
#   app <- approx(x,y,n=n)
#   x <- app$x
#   y <- app$y
#   h <- (max(x)-min(x))/(n-1)
#   if (n <= 4) stop("n needs to be >4")
#   value <- 17*(y[1]+y[n])+59*(y[2]+y[n-1])+43*(y[3]+y[n-2])+49*(y[4]+y[n-3])
#   value <- value+48*sum(y[5:(n-4)])
#   value <- (h/48)*value
#   value
# }






# main fun ----------------------------------------------------------------

#' Integrate over wavelengths
#'
#' @description Integrate over named wavelength ranges of a hyperSpec object.
#' @param x a hyperSpec object.
#' @param FUN function to be used for integration (trapezoid approximation as
#'   default).
#' @param ... (named) selection of wavelength ranges (specify as in subsetting a
#'   hyperSpec object by using a tilde).
#'
#' @return a hyperspec object with integrated spc
#' @export
#'
#' @examples
#' \dontrun{
#' integrate_wl(chondro, integral_a = 930 ~ 940, integral_b = 1230 ~ 1260)
#' }
integrate_wl <- function(x, ..., FUN = integrate_trapezoid) {
  aggregates <- list(...)
  FUN <- match.fun(FUN)
  # replace lapply with vapply here?
  tmp <- lapply(aggregates, function(.x) {
    hyperSpec::apply(
      x[, , .x], 1,
      function(.y) FUN(hyperSpec::wl(x[, , .x]), .y)
    )
  })
  out <- do.call(cbind, tmp)
  colnames(out$spc) <- names(aggregates)
  out
}
