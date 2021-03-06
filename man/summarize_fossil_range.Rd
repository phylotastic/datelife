% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_from_taxonomy.R
\name{summarize_fossil_range}
\alias{summarize_fossil_range}
\title{Summarize taxon age from PBDB to just a single min and max age}
\usage{
summarize_fossil_range(taxon, recent = FALSE, assume_recent_if_missing = TRUE)
}
\arguments{
\item{taxon}{The scientific name of the taxon you want the range of occurrences of}

\item{recent}{If TRUE, forces the minimum age to be zero}

\item{assume_recent_if_missing}{If TRUE, any taxon missing from pbdb is assumed to be recent}
}
\value{
a single row data.frame of max_ma and min_ma for the specimens, with rowname equal to taxon input
}
\description{
This uses the Paleobiology Database's API to gather information on the ages for all specimens of a taxon. It will also look for all descendants of the taxon. It fixes name misspellings if possible. It is basically a wrapper for get_fossil_range.
}
