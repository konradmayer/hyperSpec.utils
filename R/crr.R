
# helpers for crr ---------------------------------------------------------
## as given in the supplement of from https://doi.org/10.1016/j.chemolab.2018.06.009

#  Function to calculate modified Z Scores

modified_z_score <- function(x) {
  m <- stats::median(x, na.rm = TRUE)
  M <- stats::mad(x, na.rm = TRUE)
  z <- (x - m) / M
  z
}

# Function to annihalate spikes at locations z

fix_spikes_moving_average <- function(y, z, ma = 5) {
  n <- length(y)
  yout <- y
  z[1] <- z[n] <- 1
  spikes <- which(z == 1)
  for (i in spikes) {
    w <- seq(max(1, i - ma), min(n, i + ma))
    w <- w[z[w] == 0]
    yout[i] <- mean(y[w])
  }
  yout
}


#' Cosmic ray removal
#'
#' @description Cosmic spikes are often appearing in Raman spectra which were
#'   aquired using CCD detectors and need to be removed before further analysis.
#'   This function is implemented using the R code provided with Whitaker et al. 2018.
#' @param x an object of class hyperSpec.
#' @param threshold numeric.
#' @references Whitaker, Darren A., Hayes Kevin,
#'   A simple algorithm for despiking Raman spectra,
#'   Chemometrics and Intelligent Laboratory Systems,
#'   Volume 179, 2018, Pages 82-84, https://doi.org/10.1016/j.chemolab.2018.06.009.
#' @return an object of class hyperSpec with cosmic rays removed.
#' @export
#'
#' @examples
#' \dontrun{
#' chondro_spike <- chondro
#' # add spikes
#' chondro_spike[[20, , 800]] <- 5000
#' chondro_spike[[50, , 1400]] <- 4000
#' chondro_spike[[50, , 1200]] <- 5000
#' plot(chondro_spike)
#' # despike
#' despiked <- crr(chondro_spike)
#' plot(despiked)
#' }
crr <- function(x, threshold = 30) {

  # input validation
  if (!is_hyperSpec(x)) {
    stop("please provide a hyperspec object as x")
  }
  if (!(is.numeric(threshold) & (length(threshold) == 1L))) {
    stop("please provide a single numeric value as a threshold")
  }
  # calculate z scores
  y <- t(x[[]])
  w <- attr(x, "wavelength")
  z <- matrix(0, nrow(y) - 1, ncol(y))
  for (i in 1:ncol(y)) z[, i] <- modified_z_score(diff(y[, i]))
  z <- rbind(rep(0, ncol(y)), z)

  # detect threshol exceedence
  z_exceeded <- (abs(z) > threshold)
  spiky <- seq(ncol(z_exceeded)) [colSums(z_exceeded) > 0]

  # correct flagged spectra using the fixer function
  for (i in spiky) y[, i] <- fix_spikes_moving_average(y[, i], z_exceeded[, i])

  # apply despiked data to the hyperspec object
  x[[]] <- t(y)

  # add info on the despiked wavelengths - the approach is a bit hacky, but maybe to be changed into a more readable version in future
  # x$crr <- I(apply(z_exceeded, 2, function(.x) w[which(.x)])) # also always the wl succeeding the spike is flagged, therefore the lines below are introduced to fix that

  if ("crr" %in% names(x@data)) {
    warning(paste0("crr already existed in ", deparse(quote(x)), "@data and is overwritten with new values in the returned object"))
  }
  simplelag <- function(x) {
    c(NA, x[-length(x)])
  }
  x$crr <- apply(z_exceeded, 2, function(.x) w[which(.x)[which(which(.x) - simplelag(which(.x)) == 1) - 1]])
  x
}
