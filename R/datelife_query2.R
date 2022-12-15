#' Go from taxon names to a `datelifeQuery` object
#'
## #' @description
#'
#' @param input Taxon names as one of the following:
#' \describe{
#' 	 \item{A character vector of taxon names}{With taxon names as a single comma separated starting or concatenated with [c()].}
#' 	 \item{A phylogenetic tree with taxon names as tip labels}{As a `phylo` or `multiPhylo`
#' 	 			object, OR as a newick character string.}
#' }
# #' @param approximate_match Boolean; default to `TRUE`: use a slower TNRS to
# #'   correct misspellings, increasing the chance of matches (including false matches).
#' @param get_spp_from_taxon Whether to search ages for all species belonging to a
#'   given taxon or not. Default to `FALSE`. If `TRUE`, it must have same length as input.
#'   If input is a newick string with some clades it will be converted to a `phylo`
#'   object, and the order of `get_spp_from_taxon` will match `phy$tip.label`.
#' @inheritParams tnrs_match
#' @inheritDotParams rotl::tnrs_match_names -names
#' @return A `datelifeQuery` object, which is a list of four elements:
#' \describe{
#' 	 \item{$input_names}{A character vector of input taxon names.}
#' 	 \item{$tnrs_names}{A character vector of taxon names processed with TNRS.}
#' 	 \item{$ott_ids}{A numeric vector of OTT ids.}
#' 	 \item{$phy}{A `phylo` object or `NA`, if input is not a tree.}
#' }
#' @details It processes `phylo` objects and newick character string inputs
#'   with [input_process()]. If `input` is a `multiPhylo` object, only the first `phylo`
#'   element will be used. Similarly, if an `input` newick character string has multiple trees,
#'   only the first one will be used.
#' @export
make_datelife_query2 <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                                get_spp_from_taxon = FALSE,
                                reference_taxonomy = "ott",
                                ...) {
  ##############################################################################
  # enhance: add mapped (has tnrs been performed?) and matched (was it matched successfully?) element to phylo object
  # add one for each taxonomy queried: ott, catalogue of life (also contains fossils), enciclopedia of life (common names)
  if (suppressMessages(is_datelife_query(input))) {
    message("Input is already a 'datelifeQuery' object - returning input.")
    return(input)
  }
  ##############################################################################
  # input_process determines if input is newick and transforms it to phylo
  # if input is phylo or multiphylo it will also check if it is a good tree
  # if input is anything else, it will return NA
  phy_new <- input_process(input = input)
  if (inherits(phy_new, "phylo")) { # if input IS phylo
    input <- phy_new$tip.labels
  }
  ##############################################################################
  # make the datelife query
  # first, TNRS match:
  tnrs_result <- get_tnrs_names(input = input,
                 reference_taxonomy = reference_taxonomy,
                 ...)
  # second, add phy to results of TNRS match
  datelife_query <- c(tnrs_result,
                      phy = phy_new)
  # third, add the class
  class(datelife_query) <- "datelifeQuery"
  ##############################################################################
  # when get spp from taxa
  # function datelife_query_get_spp return a datelifeQuery object already
  if (any(get_spp_from_taxon)) {
    datelife_query <- datelife_query_get_spp(datelife_query,
                         get_spp_from_taxon = get_spp_from_taxon,
                         reference_taxonomy = reference_taxonomy)
  }
  ##############################################################################
  # return
  return(datelife_query)
}
