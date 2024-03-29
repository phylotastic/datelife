% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_to_table.R
\name{matrices_to_table}
\alias{matrices_to_table}
\title{Go from a list of patristic distance matrix to a table of node ages}
\usage{
matrices_to_table(matrices)
}
\arguments{
\item{matrices}{A names list of patristic distance matrices. Names correspond to the study reference.}
}
\value{
A single \code{data.frame} of "taxonA", "taxonB", and "age".
}
\description{
Go from a list of patristic distance matrix to a table of node ages
}
