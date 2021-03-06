% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/als.R
\name{als}
\alias{als}
\title{Spectral unmixing: Alternating least squares multivariate curve resolution (MCR-ALS)}
\usage{
als(x, ncomp = NULL, basis_init = NULL, prefix = "basis")
}
\arguments{
\item{x}{a hyperSpec object.}

\item{ncomp}{integer; number of pure components, does not need to be specified
if basis_init is provided.}

\item{basis_init}{a hyperSpec object; initial estimates of pure spectral components,
if basis_init is NULL, initial components are estimated
using \code{\link{opa}}.}

\item{prefix}{character; a prefix to name the pure spectra.}
}
\value{
a list with the following components:
  \describe{
    \item{coefficients}{coefficient matrix}
    \item{basis}{a hyperSpec object containing the basis (component) spectra}
    \item{summary_stats}{a named vector containing lof, rms and r2}
    \item{fit}{the als fit as returned by the custom wrapper around `ALS::als()`}
  }
}
\description{
This is a wrapper around a modified version of \link[alsace]{doALS}
  from the package {alsace} to be used with objects of class hyperSpec.
}
\seealso{
\link[alsace]{doALS}, \link[ALS]{als}
}
