% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datelife_sdm.R
\name{make_sdm}
\alias{make_sdm}
\title{Make a Super Distance Matrix (SDM) from a list of good matrices obtained with \code{\link[=get_goodmatrices]{get_goodmatrices()}}}
\usage{
make_sdm(unpadded.matrices, weighting = "flat")
}
\arguments{
\item{unpadded.matrices}{A list of patristic matrices, a \code{datelifeResult} object.}

\item{weighting}{A character vector indicating how much weight to give to each
tree in \code{input} during the SDM analysis. Options are:
\describe{
\item{weighting = "flat"}{All trees have equal weighting.}
\item{weighting = "taxa"}{Weight is proportional to number of taxa.}
\item{weighting = "inverse"}{Weight is proportional to 1 / number of taxa.}
}
Defaults to \code{weighting = "flat"}.}
}
\value{
A matrix.
}
\description{
Make a Super Distance Matrix (SDM) from a list of good matrices obtained with \code{\link[=get_goodmatrices]{get_goodmatrices()}}
}
