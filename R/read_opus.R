
# helpers -----------------------------------------------------------------

between <- function(x, left, right) {
  (x >= left) & (x <= right)
}

return_between <- function(x, left, right) {
  x[between(x, left, right)]
}

lag1 <- function(x, pad = NA) {
  stopifnot(length(pad) == 1)
  c(pad, x[-length(x)])
}

lead1 <- function(x, pad = NA) {
  stopifnot(length(pad) == 1)
  c(x[-1], pad)
}

opus_read_infoblock <- function(file_path, pr) {
  infoblock_start <- grepRaw("T01.....(TEXT)|(STRING):", pr, all = FALSE)

  if (length(infoblock_start) == 0) {
    stop("Your Bruker OPUS file may not contain any text information data block. Try to import with the argument 'read_info' set to FALSE.")
  }

  infoblock_end <- grepRaw("END", pr, offset = infoblock_start)


  stringnottext <- length(grepRaw("STRING", pr, all = TRUE)) > length(grepRaw("TEXT", pr, all = TRUE))
  infonames_start <- return_between(
    grepRaw("T\\d{2}.....(TEXT)|(STRING):",
      pr,
      all = TRUE
    ),
    infoblock_start, infoblock_end
  ) + ifelse(stringnottext, 6, 12)
  infonames_end <- vapply(
    infonames_start, function(o) {
      grepRaw("(I\\d{2})|(END)|(T\\d{2})", pr, all = FALSE, offset = o)
    },
    numeric(1)
  ) + 7


  nbytes_infonames <- infonames_end - infonames_start


  INFO_names <- vapply(seq_along(infonames_start), function(i) {
    hexView::blockString(hexView::readRaw(file_path,
      offset = infonames_start[i],
      nbytes = nbytes_infonames[i],
      human = "char", size = 1,
      endian = "little"
    ))
  }, character(1))




  info_start <- vapply(
    infonames_start, function(o) {
      tmp <- grepRaw("I\\d{2}", pr,
        all = FALSE,
        offset = o
      )
      ifelse(length(tmp) > 0, tmp, infoblock_end)
    },
    numeric(1)
  ) + 7



  fields_empty <- info_start == lead1(info_start, 0)
  if (utils::tail(info_start, 1) == (infoblock_end + 7)) {
    fields_empty[length(fields_empty)] <- TRUE
  }

  info_start <- info_start[!fields_empty]
  info_end <- vapply(
    info_start, function(o) {
      grepRaw("\\.*((T\\d{2})|(END)|(HIS))",
        pr,
        all = FALSE,
        offset = o
      )
    },
    numeric(1)
  )
  nbytes_info <- nbytes_info <- info_end - info_start

  INFO <- stats::setNames(vector("list", length(INFO_names)), nm = INFO_names)
  INFO[] <- NA

  INFO[!fields_empty] <- vapply(seq_along(info_start), function(i) {
    hexView::blockString(hexView::readRaw(file_path,
      offset = info_start[i],
      nbytes = nbytes_info[i],
      human = "char", size = 1,
      endian = "little"
    ))
  }, character(1))
  data.frame(INFO, stringsAsFactors = FALSE)
}



opus_read_spectrum <- function(file_path, pr, scale_y = TRUE) {


  # each OPUS file contains a header block in the beginning of the file

  # the integer at offset 20 specifies the number of data blocks contained
  # (including the header block)
  n_blocks_header <- hexView::readRaw(file_path,
    offset = 5 * 4, nbytes = 4,
    human = "int", size = 4
  )[[5]]

  # after the number of data blocks, for each data block 3 values are given, with
  # the third giving the start offset of the block
  header_positions <- 4 * (5 + seq_len(n_blocks_header) * 3)
  block_starts <- vapply(
    header_positions, function(p) {
      hexView::readRaw(file_path,
        offset = p, nbytes = 4,
        human = "int", size = 4
      )[[5]]
    },
    integer(1)
  )


  # read spectrum name
  snm <- grepRaw("SNM", pr, all = TRUE)[1] + 7
  snm_end <- grepRaw("END", pr, offset = snm, all = FALSE) + 11
  if (!is.na(snm)) {
    nbytes_snm <- snm_end - snm
    SNM <- hexView::blockString(hexView::readRaw(file_path,
      offset = snm,
      nbytes = nbytes_snm,
      human = "char", size = 1,
      endian = "little"
    ))
  } else {
    SNM <- NA
    warning(paste(file_path, "did not contain a name for the spectrum"))
  }

  #############################
  ### read spectrum
  ############################

  # start positions
  npt_start <- grepRaw("NPT", pr, all = TRUE) + 3 # number of discrete points
  fxv_start <- grepRaw("FXV", pr, all = TRUE) + 7 # wavelength start
  lxv_start <- grepRaw("LXV", pr, all = TRUE) + 7 # wavelength end

  # read scaling factor
  if (scale_y) {
    csf_start <- grepRaw("CSF", pr, all = TRUE) + 7 # y scaling
    CSF <- hexView::readRaw(file_path,
      width = NULL, offset = csf_start,
      nbytes = 8, human = "real", size = 8,
      endian = "little"
    )[[5]][1]
  }


  # read number of points of the spectrum
  NPT <- hexView::readRaw(file_path,
    offset = npt_start, nbytes = 12,
    human = "int", size = 4
  )[[5]][2]

  # find the spectrum data block as the one big enough to contain NPT*4 bytes
  spc_start <- block_starts[diff(block_starts) >= 4 * min(NPT)]
  spc <- hexView::readRaw(file_path,
    width = NULL, offset = spc_start,
    nbytes = NPT * 4, human = "real", size = 4,
    endian = "little"
  )[[5]]


  if (scale_y) {
    spc <- spc * CSF
  }

  # start of wavelength axis
  FXV <- hexView::readRaw(file_path,
    offset = fxv_start,
    nbytes = 16, human = "real", size = 8
  )[[5]][1]

  # end of wavelength axis
  LXV <- hexView::readRaw(file_path,
    offset = lxv_start,
    nbytes = 16, human = "real", size = 8
  )[[5]][1]

  # compute evenly spaced wavelength axis
  wavenumbers <- rev(seq(
    LXV, FXV,
    length.out = NPT
  ))

  list(wavenumbers = wavenumbers, spectrum = spc, spcname = SNM)
}


# opus_read_structure <- function(file_path, pr) {
# # TODO
# }

# main function -----------------------------------------------------------

#' Import a single Bruker OPUS binary file
#'
#' @description Function to import single OPUS binary files; supports
#'   reading meta data stored in the text information data block. At the current
#'   state, no other information (e.g. instrumental parameters) such as INS,
#'   LWN, HUM, SRC, BMS, ZFF, DAT, TIM are processed. Consider this function as
#'   experimental.
#' @param file_path character; path to the Bruker OPUS binary file.
#' @param scale_y logical; use the y scaling factor provided within the
#'   binary file.
#' @param read_info logical; read info stored in the text information data
#'   block.
#'
#' @return a hyperSpec object.
#' @export

read_opus_single <- function(file_path, scale_y = TRUE, read_info = FALSE) {

  # input validation
  if (!file.exists(file_path)) {
    stop("File does not exists. Please provide a valid path")
  }

  if (length(file_path) != 1) {
    stop("Please provide a single file path.")
  }

  # main function
  pr <- hexView::readRaw(file_path,
    offset = 0,
    nbytes = file.info(file_path)$size,
    human = "char", size = 1, endian = "little"
  )$fileRaw

  # read spectrum
  spc <- opus_read_spectrum(file_path, pr, scale_y)

  ## read info block
  if (read_info) {
    info <- opus_read_infoblock(file_path, pr)
    df <- data.frame(spcname = spc$spcname, info, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(spcname = spc$spcname, stringsAsFactors = FALSE)
  }

  suppressWarnings(
    methods::new("hyperSpec",
      spc = t(matrix(spc$spectrum)),
      wavelength = spc$wavenumbers,
      data = df
    )
  )
}



#' Read Bruker OPUS binary files
#'
#' @description This function is a wrapper around \code{\link{read_opus_single}}
#'   to import multiple OPUS binary files; supports
#'   reading meta data stored in the text information data block. At the current
#'   state, no other information (e.g. instrumental parameters) such as INS,
#'   LWN, HUM, SRC, BMS, ZFF, DAT, TIM are processed. Consider this function as
#'   experimental and sanity check imported spectra.
#' @inheritParams read_opus_single
#' @details When files to import don't share common wavelength vectors, a
#'   collapsed hyperSpec object (having a wavelength vector merged from all the
#'   individual vectors) may become very large with a lot of missing
#'   values. \code{read_opus()} can interpolate the individual spectra onto a
#'   common wavelength vector prior collapsing, either by choosing the
#'   wavelength vector of one of the imported files (giving a single integer) or
#'   by providing a new wavelength vector. Alternatively, if \code{wl} is
#'   \code{NULL}, the merged wavelength vector is used, potentially resulting in
#'   a sparse spectral matrix. If set to "interpolate" the missing values are
#'   interpolated after using the merged vector.
#' @return a hyperSpec object.
#' @param file_paths a character vector of file paths to Bruker OPUS binary
#'   files
#' @param wl numeric; either of length one giving the index of the spectrum to
#'   use the associated wavelength vector, a numeric vector of length >1
#'   specifying a new wavelength vector, \code{NULL} to use a merged wavelength
#'   vector or "interpolate" to interpolate missing values after using the
#'   merged vector.
#' @aliases read.opus
#' @name read_opus
#' @export
read_opus <- read.opus <- function(file_paths, scale_y = TRUE,
                                   read_info = FALSE,
                                   wl = 1) {
  tmp <- vector("list", length(file_paths))
  for (i in seq_along(file_paths)) {
    tryCatch(
      {
        tmp[[i]] <- read_opus_single(file_paths[i], scale_y, read_info)
      },
      error = function(err) {}
    )
  }

  errors <- vapply(tmp, is.null, logical(1))
  if (any(errors)) {
    warning(paste("It was not possible to read the file(s):",
      file_paths[errors],
      sep = "\n"
    ))
  }

  if ((!is.null(wl)) && (length(wl) == 1) && is.numeric(wl)) {
    if (!(wl %in% seq_along(file_paths))) {
      stop("please provide an index in the range of the spectra to load or a new wavelength vector")
    }
    tmp[!errors] <- lapply(
      tmp[!errors],
      function(.x) {
        hyperSpec::spc.smooth.spline(.x,
          newx = hyperSpec::wl(tmp[[wl]]),
          all.knots = TRUE
        )
      }
    )
  } else if ((!is.null(wl)) & (length(wl) > 1)) {
    tmp[!errors] <- lapply(
      tmp[!errors],
      function(.x) {
        hyperSpec::spc.smooth.spline(.x,
          newx = wl,
          all.knots = TRUE
        )
      }
    )
  }

  out <- hyperSpec::collapse(tmp[lengths(tmp) > 0L], collapse.equal = FALSE)

  if ((!is.null(wl)) && (wl == "interpolate")) {
    out <- hyperSpec::orderwl(out)
    out <- hyperSpec::spc.NA.approx(out)
  }

  out
}
