#' Get all species belonging to a taxon from the Open Tree of Life Taxonomy
#'
#' @param taxon_name A character vector providing an inclusive taxonomic name.
#' @param ott_id A numeric vector providig an Open Tree Taxonomic id number for a taxonomic name. If provided, `taxon_name` is ignored.
#'    used by Open Tree of Life Taxonomy to detect invalid taxon names.
#' @param synth_tree_only Whether to include species that are in synthetic Open Tee only or not. Default to `TRUE`.
#' @return A list of unique OTT names and OTT ids of species withing the provided taxon.
#' @export
get_opentree_species <- function(taxon_name, ott_id, synth_tree_only = TRUE) {
  if (missing(ott_id)) {
    if(missing(taxon_name)) {
      stop("A `taxon name` or an `ott id` must be provided.")
    }
    taxon_tnrs <- rotl::tnrs_match_names(names = taxon_name)
    ott_id <- taxon_tnrs$ott_id
    taxon_name <- taxon_tnrs$unique_name
    message("Getting species for OTT taxon '", 
            taxon_name,
            "' with OTT id number = ",
            ott_id, ".")
  }
  children_names <- rotl::taxonomy_subtree(ott_id = ott_id,
                               output_format = "taxa",
                               label_format = "name")$tip_label
  children_ott_ids <- rotl::taxonomy_subtree(ott_id = ott_id,
                               output_format = "taxa",
                               label_format = "id")
  children_ott_ids_numeric <- as.numeric(gsub("ott", "", children_ott_ids$tip_label))
  message("... Obtaining taxonomic info for all children within this taxon might take a while. Appreciate your patience.")
  taxon_info <- rotl::taxonomy_taxon_info(ott_ids = children_ott_ids_numeric)
  # Get unique ranks:
  message("... OTT children for this taxon belong to taxonomic ranks: ", paste0(unique(unlist(sapply(taxon_info, "[", "rank"))), collapse = ", "), ".")
  # Get subspecies
  are_subspecies <- unlist(sapply(taxon_info, "[", "rank")) == "subspecies"
  message("... There are ", sum(are_subspecies), " subspecies.")
  if(sum(are_subspecies) > 0) {
     subspecies <- unname(unlist(sapply(taxon_info[are_subspecies], "[", "name")))
    # Get unique species names from subspecies list:
    split <- strsplit(subspecies, split = " ") # split character string by blank spaces
    binomial <- lapply(split, "[", -3) # remove the subspecies epithet
    species_from_subspecies <- unique(sapply(binomial, paste0, collapse = " "))
    length(species_from_subspecies)
    # Get ott_ids of species from subspecies
    spp_from_sub_ott_ids <- rotl::tnrs_match_names(species_from_subspecies)$ott_id
    # Get their taxon info
       message("... Obtaining taxonomic info for missing species within the taxon might take a while. Appreciate your patience.")
    spp_from_sub_taxon_info <- rotl::taxonomy_taxon_info(spp_from_sub_ott_ids)
    #Unify all taxon_info
    taxon_info <- c(taxon_info, spp_from_sub_taxon_info)
  }
  ott_ids_all <-  unname(unlist(sapply(taxon_info, "[", "ott_id")))
  is_in_synth_all <- rotl::is_in_tree(ott_ids = ott_ids_all)
  are_species_all <- unlist(sapply(taxon_info, "[", "rank")) == "species"
  message("... There are ", sum(are_species_all), " species.")
  spp_in_synth <- is_in_synth_all & are_species_all
  message("... There are ", sum(spp_in_synth), " species of '", taxon_name, "' in the OpenTree synthetic tree.")
  return_species <- unlist(ifelse(synth_tree_only,
                           list(spp_in_synth),
                           list(are_species_all)))
  species_names_all <- unlist(sapply(taxon_info[return_species], "[", "name"))
  return(list(species_names = species_names_all, ott_ids = ott_ids_all))
}
