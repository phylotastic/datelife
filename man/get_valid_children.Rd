% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opentree_taxonomy.R
\name{get_valid_children}
\alias{get_valid_children}
\title{Extract valid children from given taxonomic name(s) or Open Tree of Life
Taxonomic identifiers (OTT ids) from a taxonomic source.}
\usage{
get_valid_children(input = NULL, ott_ids = NULL, reference_taxonomy = "ncbi")
}
\arguments{
\item{input}{Optional. A character vector of names or a \code{datelifeQuery} object.}

\item{ott_ids}{If not NULL, it takes this argument and ignores input. A
numeric vector of ott ids obtained with \code{\link[rotl:taxonomy_taxon_info]{rotl::taxonomy_taxon_info()}} or
\code{\link[rotl:tnrs_match_names]{rotl::tnrs_match_names()}} or \code{\link[=tnrs_match]{tnrs_match()}}.}

\item{reference_taxonomy}{A character vector with the desired taxonomic sources.
Options are "ncbi", "gbif" or "irmng". Any other value will retrieve data
from all taxonomic sources. The function defaults to "ncbi".}
}
\value{
A named list containing valid taxonomic children of given taxonomic name(s).
}
\description{
Extract valid children from given taxonomic name(s) or Open Tree of Life
Taxonomic identifiers (OTT ids) from a taxonomic source.
}
\details{
GBIF and other taxonomies contain deprecated taxa that are not marked as such
in the Open Tree of Life Taxonomy.
We are relying mainly in the NCBI taxonomy for now.
}
\examples{
# genus Dictyophyllidites with ott id = 6003921 has only extinct children
# in cases like this the same name will be returned

tti <- rotl::taxonomy_taxon_info(6003921, include_children = TRUE)
gvc <- get_valid_children(ott_ids = 6003921)

# More examples:

get_valid_children(ott_ids = 769681) # Psilotopsida
get_valid_children(ott_ids = 56601) # Marchantiophyta
}
