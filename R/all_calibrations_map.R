#' Map calibrations to nodes of a tree topology
#'
#' @param phy A phylo object to map calibrations to its nodes
#' @param calibrations A list of calibrations extracted with get_all_calibrations function
#' @return A list with two elements
#' \describe{
#'	\item{phy}{A phylo object containing the list of calibration distributions}
#'	\item{calibrations}{A data frame}
#'	}
#' @details If input is a phylo object, it is used as backbone. If it is a character vector of taxon names, an induced OToL tree is used as backbone.
#' @export
# calibrations <- get_all_calibrations(cetaceae_phyloall)
# phy <- cetaceae_phyloall[[2]]
map_all_calibrations <- function(phy, calibrations){
    # get the coincident node numbers:
    # ape::is.binary(target_tree)
    phy$tip.label <- gsub(' ', '_', phy$tip.label) #underscores vs spaces: the battle will never end.
    calibrations$taxonA <- gsub(' ', '_', calibrations$taxonA)
    calibrations$taxonB <- gsub(' ', '_', calibrations$taxonB)
    calibrations <- calibrations[which(calibrations$taxonA %in% phy$tip.label),]
    calibrations <- calibrations[which(calibrations$taxonB %in% phy$tip.label),]
    age_column <- grep("age", tolower(names(calibrations))) # get columns with "age" in their names
    target_tree_nodes <- sapply(seq(nrow(calibrations)), function(i)
      phytools::findMRCA(tree = phy, type = "node", tips =
      as.character(calibrations[i,c("taxonA", "taxonB")])))
    target_tree_nodes <- target_tree_nodes - ape::Ntip(phy)
  	all_nodes <- sort(unique(target_tree_nodes))
  	# get the node age distribution (ages taken from the first column with "age" in its name):
  	all_ages <- lapply(all_nodes, function(i) calibrations[target_tree_nodes == i, age_column[1]])
  	# any(sapply(all_ages, is.null)) # if FALSE, all nodes have at least one calibration.
  	calibrations2 <- data.frame(MRCA = paste0("n", all_nodes), MinAge = sapply(all_ages, min),
        MaxAge = sapply(all_ages, max))
    if(all(all_nodes < ape::Ntip(phy))){
        all_nodes_numbers <- all_nodes + ape::Ntip(phy)
		node_index <- "consecutive"
	} else {
        all_nodes_numbers <- all_nodes
		node_index <- "node_number"
	}
	phy$node.label <- NULL # make sure its null, so we can rename all nodes of interest to match our labels
	phy <- tree_add_nodelabels(tree = phy, node_index = node_index)  # all nodes need to be named so make_bladj_tree runs properly
    phy$calibration_distribution <- stats::setNames(all_ages, all_nodes_numbers)
    return(list(phy = phy, calibrations = calibrations2))
}
