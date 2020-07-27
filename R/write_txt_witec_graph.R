
# helper for write.txt.Witec.Graph ---------------------------------------------
write.txt.Witec.Graph_header <- function(x, path = NULL,
                                         type = "map",
                                         graphname = NULL,
                                         scan_params = NULL,
                                         connection = NULL) {
  if (!is_hyperSpecMap(x)) {
    stop("currently only export for spatially resolved hyperSpec objects is supported")
  }
  dims <- spcmap_dim(x)
  FileName <- path
  if (is.null(graphname)) {
    GraphName <- ifelse(length(unique(x$spcname)) == 1, unique(x$spcname), stop("please provide a GraphName"))
  } else {
    GraphName <- graphname
  }
  SizeX <- dims[["x"]]
  SizeY <- dims[["y"]]
  SizeGraph <- dims[["nwl"]]
  if (!is.null(x@label$.wavelength)) {
    XAxisUnit <-
      ifelse(is.expression(x@label$.wavelength),
        as.character(x@label$.wavelength),
        x@label$.wavelength
      )
  }
  if (!is.null(x@label$spc)) {
    DataUnit <-
      ifelse(is.expression(x@label$spc),
        as.character(x@label$spc),
        x@label$spc
      )
  }

  if (!is.null(scan_params)) {
    allowed <- c(
      "ScanUnit", "ScanWidth", "ScanHeight",
      "ScanOriginX", "ScanOriginY", "ScanOriginZ"
    )
    incorrect_scan_param <- !(names(scan_params) %in% allowed)
    if (any(incorrect_scan_param)) {
      warning(paste0(
        paste(names(scan_params)[incorrect_scan_param], collapse = ", "),
        " in scan_params incorrectly named.\n",
        "Only the following entries are allowed for scan_params: ", paste(allowed, collapse = ", ")
      ))
    }
    scan_params <- scan_params[names(scan_params) %in% allowed]
  }


  sink(file = connection)
  # header
  cat("//Exported ASCII-File\n[Header]\n")
  # data
  if (!is.null(FileName)) cat("FileName = ", FileName, "\n", sep = "")
  cat("GraphName = ", GraphName, "\n", sep = "")
  cat("SizeX = ", SizeX, "\n", sep = "")
  cat("SizeY = ", SizeY, "\n", sep = "")
  cat("SizeGraph = ", SizeGraph, "\n", sep = "")
  if (!is.null(scan_params)) {
    nm_scan_params <- names(scan_params)
    for (i in seq_along(scan_params)) {
      cat(nm_scan_params[[i]], " = ", scan_params[[i]], "\n", sep = "")
    }
  }
  if (exists("XAxisUnit")) cat("XAxisUnit = ", XAxisUnit, "\n", sep = "")
  if (exists("DataUnit")) cat("DataUnit = ", DataUnit, "\n", sep = "")

  sink()
}



#' Export to Witec ASCII Graph format
#'
#' @description Export Raman Spectra/Maps in the Witec ASCII Graph format to be
#'   used by external software (e.g. ImageLab).
#' @param x a \code{hyperSpec} object
#' @param file path to the export file. Specify one template filepath, the three text files for the Witec ASCII Graph format are generated automatically with the infixed " (Header)", " (X-Axis)" and " (Y-Axis)" added before the file extension.
#' @param type currently only "map" is supported.
#' @param graphname Add a custom GraphName; If not specified the function generates it from \code{x@data$spcname} in case it is unique for all spectra.
#' @param scan_params list; contains information on the scale and origin of the scan and may contain the members "ScanUnit", "ScanWidth", "ScanHeight", "ScanOriginX", "ScanOriginY", "ScanOriginZ" (take a look at the examples in case your hyperSpec object contains spatial information in "real world units" and you want to generate this list from the data).
#' @param encoding character; specifies the encoding of non-ASCII characters (e.g. for the unit), "latin1" is correct for most external software running under Windows.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(hyperSpec)
#' library(hyperSpec.utils)
#'
#' # This example uses the chondro dataset contained in {hyperSpec}
#'
#' # prepare a list for meta info in real world units calculated from values of
#' # the hyperSpec object
#' scan_params <- list(
#'   ScanUnit = "µm",
#'   ScanWidth = diff(range(chondro$x)) + 1,
#'   ScanHeight = diff(range(chondro$y)) + 1,
#'   ScanOriginX = min(chondro$x),
#'   ScanOriginY = min(chondro$y),
#'   ScanOriginZ = 0
#' )
#'
#' # write asci graph for chondro
#' write.txt.Witec.Graph(chondro,
#'   file = "graphascii_chondro.txt",
#'   graphname = "chondro", scan_params = scan_params
#' )
#'
#' # you can now use the data in external software (e.g. ImageLab),
#' # but you can also reimport it again - note that the import function from
#' # {hyperSpec} reads in using pixel units rather than "real world units"
#' # (um in this case)
#' chondro_reimport <- read.txt.Witec.Graph("graphascii_chondro (Header).txt",
#'   encoding = "latin1", type = "map"
#' )
#' # y axis is flipped after reimport - you may correct this
#' chondro_reimport$y <- -chondro_reimport$y
#' # compare original and reimport - note the differing scales
#' # (um in the original, px in the reimport;
#' # in this case 1um is represented by 1px)
#' plotmap(chondro, xlab = "", ylab = "")
#' # real world units used here with
#' # -11.55, -4.77 as the origin and an image size of 35x25um
#' plotmap(chondro_reimport, xlab = "", ylab = "")
#' # not that pixel units are used here with pixel 1,1 being the
#' # origin and an image size of 35x25px
#' }
write.txt.Witec.Graph <- function(x, file = NULL,
                                  type = "map", graphname = NULL,
                                  scan_params = NULL, encoding = "latin1") {
  if (is.null(file)) {
    stop("filename or connection needed")
  }
  type <- match.arg(type)

  path_without_ext <- tools::file_path_sans_ext(file)
  path_ext <- tools::file_ext(file)
  header_path <- paste0(path_without_ext, " (Header).", path_ext)
  xaxis_path <- paste0(path_without_ext, " (X-Axis).", path_ext)
  yaxis_path <- paste0(path_without_ext, " (Y-Axis).", path_ext)


  # write header file
  con <- file(header_path, encoding = encoding)
  write.txt.Witec.Graph_header(x,
    path = file, type = type,
    graphname = graphname, scan_params = scan_params,
    connection = con
  )
  close(con)

  # write x axis
  con <- file(xaxis_path, encoding = encoding)
  writeLines(as.character(x@wavelength), con = con)
  close(con)

  # write y axis
  con <- file(yaxis_path, encoding = encoding)
  writeLines(as.character(as.vector(t(x[[]]))), con = con)
  close(con)
}
