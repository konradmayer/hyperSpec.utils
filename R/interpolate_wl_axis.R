
# from hyperspec manual

approxfun <- function (y, wl, new.wl){
  approx(wl, y, new.wl, method = "constant",
          ties = function(x) mean(x, na.rm = TRUE)
  )$y
}

interpolate_wl_axis <- function(x) {
  hyperSpec::apply(x, 1, approxfun,
                    wl = wl(x), new.wl = unique(wl(x)),
                    new.wavelength = "new.wl")
}
