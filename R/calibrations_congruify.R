#' Use congruification to match nodes from chronograms to node in a tree topology.
#'
#' @description \code{congruify_calibrations} get nodes of a tree topology given in
#'   `phy` that correspond to the most recent common ancestor (mrca) of taxon
#'   pairs given in `calibrations`. It uses [phytools::findMRCA()] to get mrca nodes.
#'
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param chronograms A `phylo` or `multiPhylo` object, output of [datelife_search()].
#' @return A list of two elements:
#' \describe{
#' 	\item{matched_phy}{A `phylo` object with nodes renamed to match results of
#'   the mrca search. Nodes are renamed using [tree_add_nodelabels()].}
#' 	\item{matched_calibrations}{A `matchedCalibrations` object, which is the input `calibrations`
#'    object with two additional columns storing results from the mrca search with
#'    [phytools::findMRCA()]: `$mrca_node_number` and `$mrca_node_name`.}
#' 	}
#' @details The function takes pairs of taxon names in a calibrations data frame,
#' and looks for them in the vector of tip labels of the tree. If both are present,
#' then it gets the node that represents the most recent
#' common ancestor (mrca) for that pair of taxa in the tree.
#' Nodes of input `phy` can be named or not. They will be renamed.
#' @export
congruify_calibrations <- function(phy, chronograms) {
  # set up

  ##############################################################################
  ##############################################################################
  # underscores vs spaces: the battle will never end.
  ##############################################################################
  ##############################################################################
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  citations <- names(chronograms)
  chronograms <- lapply(seq(chronograms), function(x) {
    chronograms[[x]]$tip.label <- gsub(" ", "_", chronograms[[x]]$tip.label)
  })
  names(chronograms) <- citations
  ##############################################################################
  ##############################################################################
  # Extract all calibrations as a data.frame
  ##############################################################################
  ##############################################################################


  ##############################################################################
  ##############################################################################
  # Congruify source chronograms to phy
  ##############################################################################
  ##############################################################################

  ##############################################################################
  ##############################################################################
  # Match nodelabels to congruent nodes
  ##############################################################################
  ##############################################################################



  # calibrations$taxonA <- gsub(" ", "_", calibrations$taxonA)
  # calibrations$taxonB <- gsub(" ", "_", calibrations$taxonB)
  #
  # calibrations$mrca_node_number <- mrca_nodes
  # calibrations$mrca_node_name <- paste0("n", mrca_nodes)
  # calibrations$mrca_node_name <- gsub("nNA", "NA", calibrations$mrca_node_name)
  #
  # # Order the data.frame by mrca_node_number, minage and maxage
  # if ("nodeAge" %in% colnames(calibrations)) {
  #   calibrations$MinAge <- calibrations$MaxAge <- calibrations$nodeAge
  # }
  # calibrations <- calibrations[order(mrca_nodes, calibrations$MinAge, calibrations$MaxAge), ]
  #
  # # Generate node names for 'phy'
  # # All nodes need to be named so that make_bladj_tree runs properly:
  # all_nodes <- sort(unique(mrca_nodes))
  #
  # if (all(all_nodes < ape::Ntip(phy))) {
  #   all_nodes_numbers <- all_nodes + ape::Ntip(phy)
  #   node_index <- "consecutive"
  # } else {
  #   all_nodes_numbers <- all_nodes
  #   node_index <- "node_number"
  # }
  # # now we can rename all nodes of interest to match our node labels
  # # first make sure node.label is null
  # phy$node.label <- NULL
  # phy <- tree_add_nodelabels(
  #   tree = phy,
  #   node_index = node_index
  # )
  #
  return(list(matched_phy = phy,
              matched_calibrations = structure(calibrations,
                                               class = c("data.frame",
                                                         "matchedCalibrations"))))
}
