#' Process a character vector of taxon names with TNRS
#'
#' @description `make_datelife_query2` always uses TNRS (Taxonomic Name Resolution Service
#'   to process input taxon names, to correct misspellings and
#'   taxonomic name variations with [tnrs_match()], a wrapper of [rotl::tnrs_match_names()]).
#'
#' @param input Taxon names as a character vector of taxon names. Two or more
#' names can be provided as a single comma separated string or concatenated with [c()].
# #' @param approximate_match Boolean; default to `TRUE`: use a slower TNRS to
# #'   correct misspellings, increasing the chance of matches (including false matches).
#' @inheritParams tnrs_match
#' @inheritDotParams rotl::tnrs_match_names -names
#' @return A `datelifeTNRS` object, which is a list of three elements:
#' \describe{
#' 	 \item{$cleaned_names}{A character vector of names provided as input.}
#' 	 \item{$tnrs_names}{A character vector of taxon names processed with TNRS.}
#' 	 \item{$ott_ids}{A numeric vector of Open Tree of Life Taxonomy (OTT) ids.}
#' }

#' @export
get_tnrs_names <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                           reference_taxonomy = "ott",
                           ...) {
  # enhance: add mapped (has tnrs been performed?) and matched (was it matched successfully?) element to phylo object
  # add one for each taxonomy queried: ott, catalogue of life (also contains fossils), enciclopedia of life (common names)
  if (suppressMessages(is_datelife_query(input))) {
    message("Input is already a 'datelifeQuery' object - returning input.")
    return(input)
  }
  ##############################################################################
  # clean input from commas and listing format
  input <- unlist(input) # in case input is given as a list
  # split elements by the commas:
  cleaned_names <- unlist(strsplit(input, ","))
  # clean split elements of lingering unneeded white spaces:
  cleaned_names <- stringr::str_trim(cleaned_names, side = "both")

  ##############################################################################
  # process names with tnrs
  ##############################################################################
  tnrs <- tnrs_match(input = cleaned_names,
                     reference_taxonomy = reference_taxonomy,
                     ...)
  ##############################################################################
  # format return
  ##############################################################################
  # Making sure that we are using same sep as input
  # tnrs output always has spaces
  # so only modify it to underscores if input also has underscores
  if (any(grepl(pattern = "_", x = cleaned_names))) {
      tnrs$unique_name <- gsub(pattern = " ",
                          replacement = "_",
                          x = tnrs$unique_name)
  }
  # enhance: add original_taxa vector (from get_spp_from_taxon) to output here:
  datelife_query_return <- list(cleaned_names = cleaned_names,
                                tnrs_names = tnrs$unique_name,
                                ott_ids = tnrs$ott_id)
  return(structure(datelife_query_return, class = "datelifeTNRS"))
}

datelife_query_get_spp <- function(datelife_query,
                                   get_spp_from_taxon,
                                   reference_taxonomy = "ott") {
  ##############################################################################
  # checking datelife_query argument
  ##############################################################################
  is_dl_query <- suppressMessages(is_datelife_query(input = datelife_query))
  if (!is_dl_query) {
    message("* Input is not a 'datelifeQuery' object - returning NA")
    return(NA)
  }
  ##############################################################################
  # checking get_spp_from_taxon argument
  ##############################################################################
  if (all(!get_spp_from_taxon)) {
    message("Returning input datelifeQuery.")
    return(datelife_query)
  }
  if (length(get_spp_from_taxon) == 1) {
    get_spp_from_taxon <- rep(get_spp_from_taxon, length(datelife_query$cleaned_names))
  }
  ##############################################################################
  # checking both arguments have same length
  ##############################################################################
  if (length(datelife_query$ott_ids) != length(get_spp_from_taxon)) {
    message("* Number of OTT ids in input datelifeQuery and 'get_spp_from_taxon' argument have different lengths.")
    message("Please fix and try again - returning NA.")
    return(NA)
  }
  ##############################################################################
  # ott ids (and names) that we are getting species for
  ##############################################################################
  ott_ids <- datelife_query$ott_ids[get_spp_from_taxon]
  cleaned_names <- datelife_query$cleaned_names[get_spp_from_taxon]
  ##############################################################################
  # getting species
  ##############################################################################
  message("---> Getting species for the following taxon names :",
          paste0(cleaned_names, collapse = ", "), ".")
  if ("ott" %in% reference_taxonomy) {
    species_list <- lapply(ott_ids,
                           function(x) {
                             get_opentree_species(ott_id = x,
                                                  synth_tree_only = TRUE)
                            })
    return_names <- unlist(sapply(species_list, "[", "tnrs_names"))
    return_ott_ids <- unlist(sapply(species_list, "[", "ott_ids"))
    names(return_names) <- return_ott_ids
    names(return_ott_ids) <- return_names
  } else {
    message("* Other taxonomies are not implemented yet.")
    message('Please try again setting reference_taxonomy = "ott"')
    message("Returning input datelife_query.")
    return(datelife_query)
  }
  ##############################################################################
  # get data from taxa that are not getting spp from taxon
  # rearrange names:
  if (any(!get_spp_from_taxon)) {
    message("---> Adding species data from inclusive taxon names.")
    phy <- datelife_query$phy
    cleaned_names <- datelife_query$cleaned_names
    datelife_query <- lapply(datelife_query[1:3], "[", !get_spp_from_taxon)
    datelife_query$phy <- phy
    datelife_query$cleaned_names <- c(datelife_query$cleaned_names,
                                      cleaned_names[get_spp_from_taxon])
    return_names <- c(datelife_query$tnrs_names, return_names)
    return_ott_ids <- c(datelife_query$ott_ids, return_ott_ids)
  }
  ##############################################################################
  # return object
  datelife_query$tnrs_names <- return_names
  datelife_query$ott_ids <- return_ott_ids
  return(structure(datelife_query, class = "datelifeQuery"))
}
