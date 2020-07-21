spcmap2array <- function(x) {
  stopifnot(is.hyperSpecMap(x))
  dimensions <- c(x@data[,c('x', 'y')] %>% purrr::map(n_distinct), dim(x)['nwl'])
  out <- array(as.vector(x[[]]), dim = dimensions)
  attr(out, 'wavelength') <- x@wavelength
  out
}
