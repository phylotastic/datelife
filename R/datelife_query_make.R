#' Prepare a set of given taxon names as a \code{datelifeQuery} object
#'
#'
#'
#' @param input A character vector of taxon names.
#' @param use_tnrs Boolean; default to \code{FALSE}. If \code{TRUE}, use OpenTree's services
#' to resolve names. This can dramatically improve the chance of matches, but also
#' take longer.
#' @param approximate_match Boolean; default to \code{TRUE}: use a slower TNRS to correct
#' misspellings, increasing the chance of matches (including false matches)
#' @param get_spp_from_taxon boolean vector, default to \code{FALSE}. If \code{TRUE},
#' it will get all species within \code{input} taxon names. It can be selectively
#' applied, in which case it must have the same length as \code{input} taxon names.
#' It will be applied in the same order. If \code{input} is a newick string, it will
#' is converted to a \code{phylo} object, and \code{get_spp_from_taxon} is applied
#' following the order of \code{phylo$tip.label}.
#' @export
make_datelife_query_vector <- function(input = NULL,
																			 use_tnrs = FALSE,
																			 approximate_match = TRUE,
																			 get_spp_from_taxon = FALSE) {


}
