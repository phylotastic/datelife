#' Get all species belonging to a taxon from the Open Tree of Life Taxonomy (OTT)
#'
#' @param taxon_name A character vector providing an inclusive taxonomic name.
#' @param ott_id A numeric vector providig an Open Tree Taxonomic id number for
#'    a taxonomic name. If provided, `taxon_name` is ignored.
#'    Used in the context of OTT to detect invalid taxon names.
#' @param synth_tree_only Whether to include species that are in the synthetic Open
#'    Tree of Life only or not. Default to `TRUE`.
#' @return A list of unique OTT names and OTT ids of species within the provided taxon.
#' @export
get_opentree_species <- function(taxon_name, ott_id, synth_tree_only = TRUE) {
  ##############################################################################
  # checking passed arguments
  ##############################################################################
  passed <- names(as.list(match.call())[-1])
  if (all(!c("ott_id", "taxon_name") %in% passed)) {
      message("* Both 'taxon_name' and 'ott_id' arguments are missing.")
      message("Please provide at least one - returning NA.")
      return(NA)
  }
  ##############################################################################
  # getting ott ids and tnrs names
  ##############################################################################
  if (missing(ott_id)) {
    if(length(taxon_name) > 1) {
      message("This function can only process one taxon name at a time.")
      message("Please provide just one taxon name - returning NA.")
      return(NA)
    }
    taxon_tnrs <- rotl::tnrs_match_names(names = taxon_name)
    ott_id <- taxon_tnrs$ott_id
    tnrs_name <- taxon_tnrs$unique_name
    rank <- rotl::taxonomy_taxon_info(ott_id)[[1]]$rank
    message("---> Obtaining taxonomic info for all children within OTT taxon, ",
            rank, " '",
            tnrs_name, "' (input name = ", taxon_name,
            "), with OTT id number = ",
            ott_id, ".")
  }
  ##############################################################################
  # getting tnrs name from OTT ids
  ##############################################################################
  if (missing(taxon_name)) {
    if(length(ott_id) > 1) {
      message("This function can only process one OTT id at a time.")
      message("Please provide just one OTT id number - returning NA.")
      return(NA)
    }
    info <- rotl::taxonomy_taxon_info(ott_id)[[1]]
    taxon_name <- tnrs_name <- info$unique_name
    rank <- info$rank
    message("---> Obtaining taxonomic data from lineages within ",
            rank, " '",
            tnrs_name,
            "' with OTT id number = ",
            ott_id, ".")
  }
  ##############################################################################
  # getting OTT children from ott_id
  ##############################################################################
  children_names <- rotl::taxonomy_subtree(ott_id = ott_id,
                               output_format = "taxa",
                               label_format = "name")$tip_label
  children_ott_ids <- rotl::taxonomy_subtree(ott_id = ott_id,
                               output_format = "taxa",
                               label_format = "id")
  children_ott_ids_numeric <- as.numeric(gsub("ott", "", children_ott_ids$tip_label))
  taxon_info <- rotl::taxonomy_taxon_info(ott_ids = children_ott_ids_numeric)
  ##############################################################################
  # Get unique ranks:
  message("* Found the following taxonomic ranks: ",
          paste0(unique(unlist(sapply(taxon_info, "[", "rank"))),
                 collapse = ", "),
          ".")
  ##############################################################################
  # Get info from subspecies to get their species
  # because they are (always?) excluded
  are_subspecies <- unlist(sapply(taxon_info, "[", "rank")) == "subspecies"
  message("* There are ", sum(are_subspecies), " subspecies.")
  if(sum(are_subspecies) > 0) {
    subspecies <- unname(unlist(sapply(taxon_info[are_subspecies], "[", "name")))
    # Get unique species names from subspecies list:
    split <- strsplit(subspecies, split = " ") # split character string by blank spaces
    binomial <- lapply(split, "[", 1:2) # remove the subspecies epithet
    species_from_subspecies <- unique(sapply(binomial, paste0, collapse = " "))
    # length(species_from_subspecies)
    # Get ott_ids of species from subspecies
    spp_from_sub_ott_ids <- rotl::tnrs_match_names(species_from_subspecies)$ott_id
    # Get their taxon info
    spp_from_sub_taxon_info <- rotl::taxonomy_taxon_info(spp_from_sub_ott_ids)
    #Unify all taxon_info
    taxon_info <- c(taxon_info, spp_from_sub_taxon_info)
  }
  ##############################################################################
  # logical vector testing which are species in taxon_info list
  are_species <- unlist(sapply(taxon_info, "[", "rank")) == "species"
  message("* There are ", sum(are_species), " total species.")
  ##############################################################################
  # make vectors of ott ids and names for all and for species only
  ott_ids_all <-  unname(unlist(sapply(taxon_info, "[", "ott_id")))
  ott_ids_species <- ott_ids_all[are_species]
  names_all <- unlist(sapply(taxon_info, "[", "name"))
  names_species <- names_all[are_species]
  ##############################################################################
  # test if species are in tree
  is_in_synth <- rotl::is_in_tree(ott_ids = ott_ids_species)
  message("* There are ", sum(is_in_synth), " species in OpenTree's synthetic tree.")
  if (synth_tree_only) {
    message("---> Retrieving taxonomic info for ",
            sum(is_in_synth),
            " species names within '",
            taxon_name,
            "' that are in the OpenTree synthetic tree.")
    return_names <- names_species[is_in_synth]
    return_ott_ids <- ott_ids_species[is_in_synth]
    names(return_ott_ids) <- return_names
    names(return_names) <- return_ott_ids
  } else {
    message("---> Returning taxonomic info for ",
            sum(are_species),
            " species names within ", rank, " '",
            taxon_name, "'.")
    return_names <- names_species
    return_ott_ids <- ott_ids_species
    names(return_ott_ids) <- return_names
    names(return_names) <- return_ott_ids
  }
  ##############################################################################
  return(list(tnrs_names = return_names,
              ott_ids = return_ott_ids))
}


#' Quickly get all species belonging to a taxon from the Open Tree of Life Taxonomy (OTT)
#' 
#' This is less thorough than get_open_tree_species(), but much faster. It uses the fact
#' that something has just two names (genus and species) to assume that something is a
#' single species; if it has more than two names, it is assumed to be a subspecies so
#' it goes up one level in the hierarchy. It will return the subspecies and the species.
#'
#' @param taxon_name A character vector providing an inclusive taxonomic name.
#' @param ott_id A numeric vector providig an Open Tree Taxonomic id number for
#'    a taxonomic name. If provided, `taxon_name` is ignored.
#'    Used in the context of OTT to detect invalid taxon names.
#' @return A list of unique OTT names and OTT ids of species within the provided taxon.
#' @export
get_all_descendant_species <- function(taxon_name, ott_id) {
  ##############################################################################
  # checking passed arguments
  ##############################################################################
  passed <- names(as.list(match.call())[-1])
  if (all(!c("ott_id", "taxon_name") %in% passed)) {
      message("* Both 'taxon_name' and 'ott_id' arguments are missing.")
      message("Please provide at least one - returning NA.")
      return(NA)
  }
  ##############################################################################
  # getting ott ids and tnrs names
  ##############################################################################
  if (missing(ott_id)) {
    if(length(taxon_name) > 1) {
      message("This function can only process one taxon name at a time.")
      message("Please provide just one taxon name - returning NA.")
      return(NA)
    }
    taxon_tnrs <- rotl::tnrs_match_names(names = taxon_name)
    ott_id <- taxon_tnrs$ott_id
    tnrs_name <- taxon_tnrs$unique_name
    rank <- rotl::taxonomy_taxon_info(ott_id)[[1]]$rank
    message("---> Obtaining taxonomic info for all children within OTT taxon, ",
            rank, " '",
            tnrs_name, "' (input name = ", taxon_name,
            "), with OTT id number = ",
            ott_id, ".")
  }
  ##############################################################################
  # getting tnrs name from OTT ids
  ##############################################################################
  if (missing(taxon_name)) {
    if(length(ott_id) > 1) {
      message("This function can only process one OTT id at a time.")
      message("Please provide just one OTT id number - returning NA.")
      return(NA)
    }
    info <- rotl::taxonomy_taxon_info(ott_id)[[1]]
    taxon_name <- tnrs_name <- info$unique_name
    rank <- info$rank
    message("---> Obtaining taxonomic data from lineages within ",
            rank, " '",
            tnrs_name,
            "' with OTT id number = ",
            ott_id, ".")
  }
  ##############################################################################
  # getting OTT children from ott_id
  ##############################################################################
  children_names<- rotl::taxonomy_subtree(ott_id = ott_id,
								output_format = "phylo",
								label_format = "name")$tip.label
					
  children_ott_ids <- rotl::taxonomy_subtree(ott_id = ott_id,
								output_format = "phylo",
								label_format = "id")$tip.label

  children_ott_ids_numeric <- as.numeric(gsub("ott", "", children_ott_ids))
  children_names_cleaned <- gsub('_sp._|_cf._|_aff._', '_', children_names)
  taxa_with_potential_subspecies <- which(grepl("_.+_", children_names_cleaned))
  parent_species <- c()
  for(taxon_with_potential_subspecies in children_names_cleaned[taxa_with_potential_subspecies]) {
	parent_species <- c(parent_species, paste(strsplit(taxon_with_potential_subspecies, "_")[[1]][1:2], collapse = "_"))
  }
  parent_species <- unique(parent_species)
  parent_species <- parent_species[!parent_species %in% children_names]
  ottid_from_parent_species <- rotl::tnrs_match_names(parent_species)$ott_id

  good_ones <- !is.na(ottid_from_parent_species) & grepl("_", parent_species)
  parent_species <- parent_species[good_ones]
  ottid_from_parent_species <- ottid_from_parent_species[good_ones]
  
  tnrs_names = c(children_names, parent_species)
  ott_ids = gsub("ott", "", c(children_ott_ids, ottid_from_parent_species))
  good_ones_again <- grepl("_", tnrs_names) & !base::duplicated(tnrs_names)
  tnrs_names <- tnrs_names[good_ones_again]
  ott_ids <- ott_ids[good_ones_again]
  names(tnrs_names) <- ott_ids
  names(ott_ids) <- tnrs_names


  return(list(tnrs_names=tnrs_names, ott_ids=ott_ids))
}
