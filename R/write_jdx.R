

# helper ------------------------------------------------------------------


write_jdx_single <- function(wl, intensities, connection = stdout(),
                             title = NULL, jdx_version = 4.24,
                             data_type = NULL, origin = NULL, owner = NULL,
                             x_units = "1/CM", y_units = NULL,
                             ymax = 32767, xmax = NULL,
                             blocksize = 10, additional_labels = NULL) { # ymax is adviced as 32767 in format definition - scaling to an interval of -ymax:ymax, same for x

  # make sure data starts with valid value (as requested by the format definition)
  start_idx <- min(which(!is.na(intensities)))
  idx <- start_idx:length(wl)

  wl <- wl[idx]
  intensities <- intensities[idx]

  # calculate parameters
  wl_rng <- range(wl)
  firstx <- wl_rng[1]
  lastx <- wl_rng[2]
  npts <- length(wl)
  resolution <- (lastx - firstx) / (npts - 1)

  # using  packed form (PAC) as defined in MCDONALD and WILKS 1988

  # scale abscissa
  if (is.null(xmax)) {
    xfactor <- 1
    x <- (wl * xfactor)
  } else {
    xfactor <- max(abs(wl)) / xmax # see how many digits are missing to 8 digits for an integer and create multiplier
    x <- as.integer(wl / xfactor)
  }

  # scale ordinates
  yfactor <- max(abs(intensities)) / ymax # scale y to -32767 and 32767 as adviced by the file format definiton
  y <- as.integer(intensities / yfactor)



  # header

  cat("##TITLE=", title, "\n", sep = "", file = connection, append = TRUE)
  cat("##JCAMP-DX=", jdx_version, "\n", sep = "", file = connection, append = TRUE)
  cat("##DATA TYPE=", data_type, "\n", sep = "", file = connection, append = TRUE)
  cat("##ORIGIN=", origin, "\n", sep = "", file = connection, append = TRUE)
  cat("##OWNER=", owner, "\n", sep = "", file = connection, append = TRUE)

  if (!is.null(additional_labels)) {
    nm_additional_labels <- toupper(names(additional_labels))
    additional_labels <- toupper(additional_labels)
    for (i in seq_along(additional_labels)) {
      cat("##", nm_additional_labels[i], "=", additional_labels[i], "\n", sep = "", file = connection, append = TRUE)
    }
  }

  cat("##XUNITS=", x_units, "\n", sep = "", file = connection, append = TRUE)
  cat("##YUNITS=", y_units, "\n", sep = "", file = connection, append = TRUE)
  cat("##XFACTOR=", xfactor, "\n", sep = "", file = connection, append = TRUE)
  cat("##YFACTOR=", yfactor, "\n", sep = "", file = connection, append = TRUE)
  cat("##FIRSTX=", firstx, "\n", sep = "", file = connection, append = TRUE)
  cat("##LASTX=", lastx, "\n", sep = "", file = connection, append = TRUE)
  cat("##RESOLUTION=", resolution, "\n", sep = "", file = connection, append = TRUE)
  cat("##MAXY=", ymax, "\n", sep = "", file = connection, append = TRUE)
  cat("##MINY=", -ymax, "\n", sep = "", file = connection, append = TRUE)
  cat("##NPOINTS=", npts, "\n", sep = "", file = connection, append = TRUE)
  cat("##FIRSTY=", intensities[1], "\n", sep = "", file = connection, append = TRUE)
  cat("##XYDATA=", "(X++(Y..Y))", "\n", sep = "", file = connection, append = TRUE)

  # tabular spectral data
  yprint <- formatC(y, flag = "+")
  blockends <- seq(from = blocksize, to = length(wl), by = blocksize)

  for (be in blockends) {
    cat(wl[be - blocksize + 1], yprint[(be - blocksize + 1):be],
      "\n",
      sep = "", file = connection, append = TRUE
    )
  }

  # if number of values is not dividable by blocksize - add the remaining incomplete bloc
  if ((length(wl) %% blocksize) != 0) {
    startlastblock <- (blockends[length(blockends)] + 1)
    remaining_idx <- startlastblock:length(wl)
    cat(wl[startlastblock], yprint[remaining_idx],
      "\n",
      sep = "", file = connection, append = TRUE
    )
  }

  # end block
  cat("##END=", sep = "", file = connection, append = TRUE)
}


# library(hyperSpec)
# x <- read.jdx('/media/konrad/raman2/6 Team/14 Konrad/test-readopus/xx2/532 4-Allyl-2,6-dimethoxyphenol (Avall 0+90).dx', encoding = 'latin1')
# intensities <- x$spc[1,]
# wl <- wl(x)


#' Write to file in JCAMP-DX format
#'
#' @description  Writes spectra to a single or individual file in JCAMP-DX
#'   format. The format definition for JCAMP-DX 4.24 as specified in
#'   McDonald and Wilks (1988) is used. Please refer to this publication for
#'   more information about the data labels. This implementation writes labels
#'   for the required entries (table IIa) plus missing entries for the
#'   required and optional spectral parameters (see table V).
#' @param x a hyperSpec object
#' @param path character; if only one path is provided, all spectra are written
#'   as blocks to a single file; specify a path for each spectrum to write
#'   separate files.
#' @param encoding character; usually "latin1"
#' @param title character; title for each spectrum, usually a column of the
#'   hyperSpec object is used here, such as \code{x@data$spc}.
#' @param data_type character; translated to upper case; allowed  values are
#'   NULL, "RAMAN SPECTRUM", "INFRARED SPECTRUM".
#' @param origin character; length 1 or nrow(xx)
#' @param owner character; length 1 or nrow(xx)
#' @param x_units character; translated to upper case; allowed  values are
#'   "1/CM", "MICROMETERS", "NANOMETERS".
#' @param y_units character; translated to upper case; allowed  values are
#'   "ARBITRARY UNITS", "ABSORBANCE", "TRANSMITTANCE", "REFLECTANCE".
#' @param ymax numeric; ordinate values is scaled to a range of c(-ymax, ymax)
#'   before being converted to integer.
#' @param xmax numeric; , abscissa is scaled to a range of c(-xmax, xmax) before
#'   being converted to integer. If NULL, the unconverted doubles are written
#'   to the file.
#' @param additional_labels list; provide a list of named character vectors
#'   which get used as additional header entries after being translated to
#'   uppercase characters - e.g. for two spectra
#'   \code{list(
#'   c("date" = "2020-01-01", "Sample description" = "Sample1"),
#'   c("date" = "2020-01-05", "Sample description" = "Sample2"))}. While you
#'   enter any labels here its best to stick to the ones defined in McDonald and
#'   Wilks (1988) for software compatibility.
#' @param blocksize integer; number of ordinate values in a row.
#' @references ROBERT S. McDONALD and PAUL A. WILKS, JR. 1988. JCAMP-DX: A
#'   Standard Form for the Exchange of Infrared Spectra in Computer Readable
#'   Form. Appl. Spectrosc. 42(1), pp151-162
#' @name write_jdx
#' @aliases write.jdx
#' @export

write_jdx <- write.jdx <- function(x, path = x@data$filename, encoding = "latin1",
                                   title = x@data$spcname,
                                   data_type = c(NA_character_, "RAMAN SPECTRUM", "INFRARED SPECTRUM"),
                                   origin = NA, owner = NA,
                                   x_units = c(NA_character_, "1/CM", "MICROMETERS", "NANOMETERS"),
                                   y_units = c(
                                     "ARBITRARY UNITS", "ABSORBANCE",
                                     "TRANSMITTANCE", "REFLECTANCE"
                                   ),
                                   additional_labels = NULL,
                                   ymax = 32767, xmax = NULL, blocksize = 10) {

  # input validation
  data_type <- match.arg(data_type)
  x_units <- match.arg(x_units)
  y_units <- match.arg(y_units)

  if (is.null(path)) {
    stop("please provide a valid path or a vector of paths with a length equal to the number of spectra in x.")
  }

  # repeat arguments if only lengh 1 is specified or test if length is nrow(x) if otherwise
  args <- c(
    "title", "data_type", "origin", "owner", "x_units", "y_units",
    "ymax", "xmax", "blocksize"
  )
  for (arg in args) {
    tmp <- get(arg)
    if (is.character(tmp)) {
      tmp <- toupper(tmp)
    }
    if ((length(tmp) == 1) | is.null(tmp)) {
      assign(arg, rep(tmp, nrow(x)))
    } else if (length(tmp) == nrow(x)) {
      assign(arg, tmp)
    } else {
      stop(paste("Argument", arg, "has to be of length 1 or the same length as number of spectra."))
    }
  }


  intensities_mat <- x[[]]
  wl <- hyperSpec::wl(x)

  # CASE 1: a path is provided for each spectrum
  if (length(path) > 1) {
    if (nrow(intensities_mat) != length(path)) {
      stop("provide a file path for each spectrum")
    }

    for (i in seq_len(nrow(x))) {
      con <- file(path[i], encoding = encoding)
      open(con, "w")
      write_jdx_single(wl, intensities_mat[i, ],
        connection = con,
        title = title[i],
        jdx_version = 4.24,
        data_type = data_type[i], origin = origin[i],
        owner = owner[i], x_units = x_units[i],
        y_units = y_units[i],
        ymax = ymax[i], xmax = xmax[i],
        blocksize = blocksize[i],
        additional_labels = additional_labels[[i]]
      )
      close(con)
    }
  }

  # CASE 2: only one path is provided, all spectra are written to a single file
  if (length(path) == 1) {
    con <- file(path, encoding = encoding)
    open(con, "w")
    nblocks <- seq_len(nrow(x))

    cat("##TITLE=COLLECTION OF SPECTRA", "\n", sep = "", file = con, append = TRUE)
    cat("##JCAMP-DX=", 4.24, "\n", sep = "", file = con, append = TRUE)
    cat("##DATA TYPE=LINK", "\n", sep = "", file = con, append = TRUE)
    cat("##BLOCKS=", nblocks + 1, "\n\n\n", sep = "", file = con, append = TRUE)

    for (i in nblocks) {
      write_jdx_single(wl, intensities_mat[i, ],
        connection = con,
        title = title[i],
        jdx_version = 4.24,
        data_type = data_type[i], origin = origin[i],
        owner = owner[i], x_units = x_units[i],
        y_units = y_units[i],
        ymax = ymax[i], xmax = xmax[i],
        blocksize = blocksize[i],
        additional_labels = additional_labels[[i]]
      )
      if (nrow(x) > 1) {
        cat("\n\n", sep = "", file = con, append = TRUE)
      }
    }

    cat("##END=$$ COLLECTION OF SPECTRA", "\n", sep = "", file = con, append = TRUE)
    close(con)
  }
}
