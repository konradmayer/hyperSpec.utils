is_hyperSpec <- is.hyperSpec <- function (x)
{
  inherits(x, "hyperSpec")
}



is_hyperSpecMap <- is.hyperSpecMap <- function (x)
{
  inherits(x, "hyperSpec") & all(c('x', 'y') %in% names(x@data))
}
