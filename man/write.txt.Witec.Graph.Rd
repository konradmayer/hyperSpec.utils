% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_txt_witec_graph.R
\name{write.txt.Witec.Graph}
\alias{write.txt.Witec.Graph}
\title{Export to Witec ASCII Graph format}
\usage{
write.txt.Witec.Graph(
  x,
  file = NULL,
  type = "map",
  graphname = NULL,
  scan_params = NULL,
  encoding = "latin1"
)
}
\arguments{
\item{x}{a \code{hyperSpec} object}

\item{file}{path to the export file. Specify one template filepath, the three text files for the Witec ASCII Graph format are generated automatically with the infixed " (Header)", " (X-Axis)" and " (Y-Axis)" added before the file extension.}

\item{type}{currently only "map" is supported.}

\item{graphname}{Add a custom GraphName; If not specified the function generates it from \code{x@data$spcname} in case it is unique for all spectra.}

\item{scan_params}{list; contains information on the scale and origin of the scan and may contain the members "ScanUnit", "ScanWidth", "ScanHeight", "ScanOriginX", "ScanOriginY", "ScanOriginZ" (take a look at the examples in case your hyperSpec object contains spatial information in "real world units" and you want to generate this list from the data).}

\item{encoding}{character; specifies the encoding of non-ASCII characters (e.g. for the unit), "latin1" is correct for most external software running under Windows.}
}
\description{
Export Raman Spectra/Maps in the Witec ASCII Graph format to be
  used by external software (e.g. ImageLab).
}
\examples{
\dontrun{
library(hyperSpec)
library(hyperSpec.utils)

# This example uses the chondro dataset contained in {hyperSpec}

# prepare a list for meta info in real world units calculated from values of
# the hyperSpec object
scan_params <- list(
  ScanUnit = "µm",
  ScanWidth = diff(range(chondro$x)) + 1,
  ScanHeight = diff(range(chondro$y)) + 1,
  ScanOriginX = min(chondro$x),
  ScanOriginY = min(chondro$y),
  ScanOriginZ = 0
)

# write asci graph for chondro
write.txt.Witec.Graph(chondro,
  file = "graphascii_chondro.txt",
  graphname = "chondro", scan_params = scan_params
)

# you can now use the data in external software (e.g. ImageLab),
# but you can also reimport it again - note that the import function from
# {hyperSpec} reads in using pixel units rather than "real world units"
# (um in this case)
chondro_reimport <- read.txt.Witec.Graph("graphascii_chondro (Header).txt",
  encoding = "latin1", type = "map"
)
# y axis is flipped after reimport - you may correct this
chondro_reimport$y <- -chondro_reimport$y
# compare original and reimport - note the differing scales
# (um in the original, px in the reimport;
# in this case 1um is represented by 1px)
plotmap(chondro, xlab = "", ylab = "")
# real world units used here with
# -11.55, -4.77 as the origin and an image size of 35x25um
plotmap(chondro_reimport, xlab = "", ylab = "")
# not that pixel units are used here with pixel 1,1 being the
# origin and an image size of 35x25px
}
}
