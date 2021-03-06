% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datelife_query.R
\name{input_process}
\alias{input_process}
\title{Takes a phylo object or a character string and figure out if it's correct newick format or a list of species}
\usage{
input_process(input, verbose = FALSE)
}
\arguments{
\item{input}{A character vector of taxon names, a tree as a 'phylo' object or a newick character string, or a 'datelifeQuery' object from make_datelife_query function.}

\item{verbose}{Boolean. If TRUE, it gives printed updates to the user.}
}
\value{
A phylo object or NA if no tree
}
\description{
Takes a phylo object or a character string and figure out if it's correct newick format or a list of species
}
