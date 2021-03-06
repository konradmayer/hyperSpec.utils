% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/als.R
\name{opa}
\alias{opa}
\title{Orthogonal Projection Approach}
\usage{
opa(x, ncomp = NULL, return_scaled = TRUE)
}
\arguments{
\item{x}{a hyperSpec object.}

\item{ncomp}{integer, number of most dissimilar spectra to return.}

\item{return_scaled}{if TRUE, unit-length scaled spectra are returned,
if FALSE the original intensities are retained.}
}
\value{
a hyperSpec object containing ncomp spectra.
}
\description{
find the most dissimilar spectra in a hyperSpec object. This
  function is a customized version of the function \link[alsace]{opa} from
  the package {alsace} to work with objects of the class hyperSpec.
}
