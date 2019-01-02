# functions to create a tree from taxonomy

#' Gets classification paths for a vector of taxa
#'
#' This uses the taxize package's wrapper of the Global Names Resolver to get taxonomic paths for the vector of taxa you pass in. Sources is a vector of source labels in order (though it works best if everything uses the same taxonomy, so we recommend doing just one source). You can see options by doing taxize::gnr_datasources(). Our default is Catalogue of Life
#' @param taxa Vector of taxon names
#' @param sources Vector of names of preferred sources; see taxize::gnr_datasources()
#' @return A list with resolved taxa (a tibble, from taxize::gnr_resolve) and a vector of taxa not resolved
#' @export
classification_paths_from_taxonomy <- function(taxa, sources="Catalogue of Life") {
  all_sources <- taxize::gnr_datasources()
  source_ids <- rep(NA, length(sources))
  for (i in seq_along(sources)) {
    source_ids[i] <- all_sources$id[grepl(sources[i], all_sources$title, ignore.case=TRUE)]
  }
  source_ids <- source_ids[!is.na(source_ids)]
  resolved_taxa <- taxize::gnr_resolve(taxa, best_match_only=TRUE, fields="all", preferred_data_sources=source_ids, with_context=TRUE)
  resolved_taxa <- resolved_taxa[grepl("\\|", resolved_taxa$classification_path),] # sometimes PBDB resolves taxa but doesn't have a taxonomy for them. This isn't what we want.
  missing_taxa <- taxa[!taxa %in% resolved_taxa$user_supplied_name]

  return(list(resolved=resolved_taxa, unresolved=missing_taxa))
}


#' Gets a taxonomic tree from a vector of taxa
#'
#' This uses the taxize package's wrapper of the Global Names Resolver to get taxonomic paths for the vector of taxa you pass in. Sources is a vector of source labels in order (though it works best if everything uses the same taxonomy, so we recommend doing just one source). You can see options by doing taxize::gnr_datasources(). Our default is Catalogue of Life. The output is a phylo object (typically with many singleton nodes if collapse_singles is FALSE: nodes with only one descendant (like "Homo" having "Homo sapiens" as its only descendant) but these singletons typically have node.labels
#' @param taxa Vector of taxon names
#' @param sources Vector of names of preferred sources; see taxize::gnr_datasources()
#' @param collapse_singles If true, collapses singleton nodes
#' @return A list containing a phylo object with resolved names and a vector with unresolved names
#' @export
#' @examples
#' taxa <- c("Homo sapiens", "Ursus arctos", "Pan paniscus", "Tyrannosaurus rex", "Ginkgo biloba", "Vulcan", "Klingon")
#' results <- tree_from_taxonomy(taxa)
#' print(results$unresolved) # The taxa that do not match
#' ape::plot.phylo(results$phy) # may generate warnings due to problems with singletons
#' ape::plot.phylo(ape::collapse.singles(results$phy), show.node.label=TRUE)
#' # got rid of singles, but this also removes a lot of the node.labels
tree_from_taxonomy <- function(taxa, sources="Catalogue of Life", collapse_singles=TRUE) {
  classification_results <- classification_paths_from_taxonomy(taxa=taxa, sources=sources)
  paths <- classification_results$resolved$classification_path
  if(length(paths)<2) { #not enough for a tree
    return(NULL)
  }
  paths <- gsub('Not assigned\\|', "", paths) # CoL gives Not assigned for some taxa. We don't want these. Then we split
  paths <- strsplit( paths, "\\|")
  leaves <- sapply(paths, tail, n=1)
  tip.label <- leaves
  node.label <- rev(unique(unlist(paths)))
  node.label <- node.label[!node.label %in% leaves]
  node.label <- c("Life", node.label)
  edge <- matrix(nrow=0, ncol=2)
  for(path_index in seq_along(paths)) {
    local_path <- rev(paths[[path_index]])
    tip.id <- which(local_path[1] == tip.label)
    ancestor.id <- 1+length(tip.label) + which(local_path[2] == node.label)
    edge <- rbind(edge, c(ancestor.id, tip.id))
    for(step_index in seq_along(local_path)) {
      if(step_index>1 & step_index<length(local_path)) { # since we've done the first step already, and the last node requires special treatment (since it connects to root)
        tipward.id <- 1 + length(tip.label) + which(local_path[step_index] == node.label)
        rootward.id <- 1 + length(tip.label) + which(local_path[step_index+1] == node.label)
        edge <- rbind(edge, c(rootward.id, tipward.id))
      }
      if(step_index==length(local_path)) {
        tipward.id <- 1+length(tip.label) + which(local_path[step_index] == node.label)
        edge <- rbind(edge, c(1+length(tip.label), tipward.id))
      }
    }
  }
  edge <- unique(edge) #get rid of edges added multiple times
  #edge <- edge[edge[,1]!=edge[,2],] # sometimes have the same edge be its own parent. This is bad.
  phy <- list(edge=edge, tip.label=tip.label, node.label=node.label, Nnode=nrow(edge))
  class(phy) <- "phylo"
  phy <- ape::reorder.phylo(phy)
  if(collapse_singles) {
    phy <- ape::collapse.singles(phy)
  }
  return(list(phy=phy, unresolved=classification_results$unresolved))
}

#' Get the ages for a taxon from PBDB
#'
#' This uses the Paleobiology Database's API to gather information on the ages for all specimens of a taxon. It will also look for all descendants of the taxon. It fixes name misspellings if possible.
#' @param taxon The scientific name of the taxon you want the range of occurrences of
#' @param recent If TRUE, forces the minimum age to be zero
#' @param assume_recent_if_missing If TRUE, any taxon missing from pbdb is assumed to be recent
#' @return a data.frame of max_ma and min_ma for the specimens
#' @export
get_fossil_range <- function(taxon, recent=FALSE, assume_recent_if_missing=TRUE) {
  all_sources <- taxize::gnr_datasources()
  source_id <- all_sources$id[grepl("The Paleobiology Database", all_sources$title, ignore.case=TRUE)]
  taxon_gnred_df <- taxize::gnr_resolve(taxon, best_match_only=TRUE, fields="all", preferred_data_sources=source_id)
  if(nrow(taxon_gnred_df)==0) {
    if(assume_recent_if_missing) {
      dates <- data.frame(max_ma=0, min_ma=0)
      rownames(dates) <- taxon
      return(dates)
    } else {
      return(NA)
    }
  }
  taxon_gnred <- tail(strsplit( taxon_gnred_df$classification_path, "\\|")[[1]],1)
  taxon_string <- utils::URLencode(taxon_gnred)
  dates <- data.frame()
  try(dates <- read.csv(url(paste0("https://paleobiodb.org/data1.2/occs/list.txt?base_name=", taxon_string))))
  try(dates <- dates[,c("max_ma", "min_ma")])
  if(recent) { #we know it still exists
    dates <- rbind(dates, c(min(c(0, dates$min_ma), na.rm=TRUE),0))
  } else { #check to see if exists
    taxon_info <- NULL
    try(taxon_info <- read.csv(url(paste0("https://paleobiodb.org/data1.2/taxa/single.txt?name=", taxon_string))))
    if(!is.null(taxon_info)) {
      if(taxon_info$is_extant=="extant") {
        dates <- rbind(dates, c(min(c(0, dates$min_ma), na.rm=TRUE),0)) #so that even if it has no fossils, record it exists
      }
    }
  }
  if(any(is.na(dates)) & assume_recent_if_missing) {
    dates[1,] <- c(0,0)
  }
  return(dates)
}

#' Summarize taxon age from PBDB to just a single min and max age
#'
#' This uses the Paleobiology Database's API to gather information on the ages for all specimens of a taxon. It will also look for all descendants of the taxon. It fixes name misspellings if possible. It is basically a wrapper for get_fossil_range.
#' @param taxon The scientific name of the taxon you want the range of occurrences of
#' @param recent If TRUE, forces the minimum age to be zero
#' @param assume_recent_if_missing If TRUE, any taxon missing from pbdb is assumed to be recent
#' @return a single row data.frame of max_ma and min_ma for the specimens, with rowname equal to taxon input
#' @export
summarize_fossil_range <- function(taxon, recent=FALSE, assume_recent_if_missing=TRUE) {
  dates <- get_fossil_range(taxon, recent, assume_recent_if_missing)
  result <- data.frame(max_ma=max(dates$max_ma), min_ma=min(dates$min_ma), stringsAsFactors=FALSE)
  rownames(result) <- taxon
  return(result)
}

#' Date with Paleobiology Database and paleotree
#'
#' This will take a topology, look up information about fossils for taxa on the tree, and use paleotree's timePaleoPhy() function to compute branch lengths.
#' @param phy Phylogeny of taxa
#' @param recent If TRUE, forces the minimum age to be zero for any taxon
#' @param assume_recent_if_missing If TRUE, any taxon missing from pbdb is assumed to be recent
#' @return A dated tree
#' @export
#' @examples
#' taxa <- c("Archaeopteryx", "Pinus", "Quetzalcoatlus", "Homo sapiens",
#'        "Tyrannosaurus rex", "Megatheriidae", "Metasequoia", "Aedes", "Panthera")
#' phy <- tree_from_taxonomy(taxa, sources="The Paleobiology Database")$phy
# uncomment the following part when we paleotree can be imported
# #' chronogram <- date_with_pbdb(phy)
# #' ape::plot.phylo(chronogram)
# #' ape::axisPhylo()
date_with_pbdb <- function(phy, recent=FALSE, assume_recent_if_missing=TRUE) {
  all_dates <- as.data.frame(t(sapply(phy$tip.label, summarize_fossil_range, recent=recent, assume_recent_if_missing=assume_recent_if_missing)))
  all_dates$max_ma <- as.numeric(all_dates$max_ma)
  all_dates$min_ma <- as.numeric(all_dates$min_ma)
  # enhance: paleotree cannot be included as an import, but it needs to be if we wanna use them inside our functions
  chronogram <- paleotree::timePaleoPhy(phy, all_dates, add.term=TRUE)
  return(chronogram)
}
