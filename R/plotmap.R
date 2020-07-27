#' Plotmap using viridis
#'
#' @description Use the viridis colorramp in plotmap as default. The model can
#'   contain the special column name .wavelength to specify the wavelength axis.
#'
#' @inheritParams hyperSpec::plotmap
#' @export
plotmap_viridis <- function(object, model = spc ~ x * y, func = mean,
                            func.args = list(), ...) {
  hyperSpec::plotmap(object,
    model = model, func = func, func.args = func.args,
    par.settings = rasterVis::viridisTheme(), ...
  )
}
