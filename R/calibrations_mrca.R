#' Identify nodes of a tree topology that are most recent common ancestor (mrca)
#' of taxon pairs from a `calibrations` object
#'
#' @description \code{mrca_calibrations} get nodes of a tree topology given in
#'   `phy` that correspond to the most recent common ancestor (mrca) of taxon
#'   pairs given in `calibrations`. It uses [phytools::findMRCA()] to get mrca nodes.
#'
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param calibrations A `calibrations` object, an output of
#'   [extract_calibrations_phylo()].
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
mrca_calibrations <- function(phy, calibrations) {
  # Should we implement this??? -> If input is a phylo object, it is used as backbone. If it is a character vector
  # of taxon names, an induced synthetic OpenTree subtree is used as backbone.
  # calibrations <- get_all_calibrations(cetaceae_phyloall)
  # phy <- cetaceae_phyloall[[2]]
  # get the coincident node numbers:
  # ape::is.binary(target_tree)
  if (!inherits(phy, "phylo")) {
    warning("'phy' is not a 'phylo' object.\nCalibrations can't be matched.")
    return(NA)
  }
  if (!inherits(calibrations, "data.frame")) {
    warning("'calibrations' is not a 'data.frame'.\nCalibrations can't be matched.")
    return(NA)
  }
  phy$tip.label <- gsub(" ", "_", phy$tip.label) # underscores vs spaces: the battle will never end.
  calibrations$taxonA <- gsub(" ", "_", calibrations$taxonA)
  calibrations$taxonB <- gsub(" ", "_", calibrations$taxonB)
  # identify taxon name pairs that are in phy:
  in_phy <- calibrations$taxonA %in% phy$tip.label & calibrations$taxonB %in% phy$tip.label
  if (all(!in_phy)) {
    warning("Taxon name pairs in 'calibrations' do not match any tip labels in 'phy'.")
    return(NA)
  }
  mrca_nodes <- sapply(
    seq(in_phy),
    function(i) {
      if (in_phy[i]) {
        x <- phytools::findMRCA(
                    tips = as.character(calibrations[i, c("taxonA", "taxonB")]),
                    tree = phy,
                    type = "node"
                  )
        # phytools returns node numbers starting at Ntip(phy) + 1
        # x - ape::Ntip(phy)
     } else { NA }
    }
  )

# attempt to simplify the code that gets mrca_nodes:
# tryCatch(expr = {
#   x <- phytools::findMRCA(tips = as.character(calibrations[i, c("taxonA", "taxonB")]),
#                      tree = phy,
#                      type = "node")
#   # phytools returns node numbers starting at Ntip(phy) + 1
#   x - ape::Ntip(phy)
# }, error = NA)

# taxon pairs not in 'phy':
# calibrations[which(is.na(mrca_nodes)), -6]

  calibrations$mrca_node_number <- mrca_nodes
  calibrations$mrca_node_name <- paste0("n", mrca_nodes)
  calibrations$mrca_node_name <- gsub("nNA", "NA", calibrations$mrca_node_name)

  # Order the data.frame by mrca_node_number, minage and maxage
  if ("nodeAge" %in% colnames(calibrations)) {
    calibrations$MinAge <- calibrations$MaxAge <- calibrations$nodeAge
  }
  calibrations <- calibrations[order(mrca_nodes, calibrations$MinAge, calibrations$MaxAge), ]

  # Generate node names for 'phy'
  # All nodes need to be named so that make_bladj_tree runs properly:
  all_nodes <- sort(unique(mrca_nodes))

  if (all(all_nodes < ape::Ntip(phy))) {
    all_nodes_numbers <- all_nodes + ape::Ntip(phy)
    node_index <- "from_1"
  } else {
    all_nodes_numbers <- all_nodes
    node_index <- "node_number"
  }
  # now we can rename all nodes of interest to match our node labels
  # first make sure node.label is null
  phy$node.label <- NULL
  phy <- tree_add_nodelabels(
    tree = phy,
    node_index = node_index
  )
  #
  return(list(matched_phy = phy,
              matched_calibrations = structure(calibrations,
                                               class = c("data.frame",
                                                         "mrcaCalibrations"))))
}
