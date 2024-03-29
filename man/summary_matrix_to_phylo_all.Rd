% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_phylo_all.R
\name{summary_matrix_to_phylo_all}
\alias{summary_matrix_to_phylo_all}
\title{Get minimum, median, mean, midpoint, and maximum summary chronograms from a
summary matrix of a \code{datelifeResult} object.}
\usage{
summary_matrix_to_phylo_all(
  summ_matrix,
  datelife_query = NULL,
  target_tree = NULL,
  total_distance = TRUE,
  ...
)
}
\arguments{
\item{summ_matrix}{Any summary patristic distance matrix, such as the ones obtained with \code{\link[=datelife_result_sdm_matrix]{datelife_result_sdm_matrix()}} or \code{\link[=datelife_result_median_matrix]{datelife_result_median_matrix()}}.}

\item{datelife_query}{A \code{datelifeQuery} object, usually an output of \code{\link[=make_datelife_query]{make_datelife_query()}}.}

\item{target_tree}{A \code{phylo} object. Use this in case you want a specific
backbone for the output tree.}

\item{total_distance}{Whether the input \code{summ_matrix} stores total age distance
(from tip to tip) or distance from node to tip. Default to \code{TRUE},
divides the matrix in half, if \code{FALSE} it will take it as is.}

\item{...}{
  Arguments passed on to \code{\link[=get_otol_synthetic_tree]{get_otol_synthetic_tree}}
  \describe{
    \item{\code{otol_version}}{Version of Open Tree of Life to use}
    \item{\code{resolve}}{Defaults to \code{TRUE}. Whether to resolve the tree at random or not.}
    \item{\code{input}}{Optional. A character vector of names or a \code{datelifeQuery} object.}
    \item{\code{ott_ids}}{If not NULL, it takes this argument and ignores input. A
numeric vector of ott ids obtained with \code{\link[rotl:taxonomy_taxon_info]{rotl::taxonomy_taxon_info()}} or
\code{\link[rotl:tnrs_match_names]{rotl::tnrs_match_names()}} or \code{\link[=tnrs_match]{tnrs_match()}}.}
  }}
}
\value{
A \code{multiPhylo} object of length 5. It contains min, mean, median, midpoint, and max summary chronograms.
}
\description{
Get minimum, median, mean, midpoint, and maximum summary chronograms from a
summary matrix of a \code{datelifeResult} object.
}
\details{
With this function users can choose the minimum, mean or maximum ages from
the summary matrix as calibration points to get a single summary chronogram.
Users get all three summary chronograms in a \code{multiPhylo} object.
}
