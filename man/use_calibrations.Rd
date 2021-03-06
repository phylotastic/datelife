% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_calibrations_use.R
\name{use_calibrations}
\alias{use_calibrations}
\title{Perform a dating analysis on a tree topology using a given set of calibrations.}
\usage{
use_calibrations(
  phy = NULL,
  calibrations = NULL,
  dating_method = "bladj",
  type = "median",
  ...
)
}
\arguments{
\item{phy}{A phylo object with or wothput branch lengths}

\item{calibrations}{A data frame of calibrations from get_all_calibrations function, or a subset of it.}

\item{dating_method}{The method used for tree dating.}

\item{type}{The type of age to use as calibration: "median", "mean", "min", or "max".}

\item{...}{
  Arguments passed on to \code{\link[=use_calibrations_pathd8]{use_calibrations_pathd8}}
  \describe{
    \item{\code{expand}}{How much to expand by each step to get consistent calibrations. Should be between 0 and 1.}
    \item{\code{giveup}}{How many expansions to try before giving up}
  }}
}
\value{
A phylo object with branch lengths proportional to time.
}
\description{
Perform a dating analysis on a tree topology using a given set of calibrations.
}
\details{
If phy does not have branch lengths, dating_method is ignored and BLADJ will be used.
}
