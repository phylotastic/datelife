\name{nonmatching_by_row}
\alias{nonmatching_by_row}
\title{Find only rows in a data.frame where the values don't match in a row}
\usage{
  nonmatching_by_row(treename, cleanedtree)
}
\description{
  Find only rows in a data.frame where the values don't
  match in a row
}
\examples{
dff <- data.frame(FabreEtAl2009[[1]][[1]]$tip.label, out$tip.label)
nonmatching_by_row(treename=dff)
}

