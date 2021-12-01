#' Match calibrations to nodes of a tree
#'
#' @description \code{match_all_calibrations} summarizes nodes in tree that match
#' to any pair of given taxon names in a data frame of secondary calibrations.
#'
#' @param phy A \code{phylo} object. Nodes can be named or not.
# #' or a vector of taxon names (see details).
#' @param calibrations An object of class \code{datelifeCalibrations}, that is,
#' an output of [get_all_calibrations()].
#' @return A list with:
#' \describe{
#' 	\item{phy}{A phylo object with additional data of calibration distributions per node.}
#' 	\item{matched_calibrations}{A data frame of summarized calibrations.}
#' 	\item{present_calibrations}{A data frame of summarized calibrations.}
#' 	}
#' TODO: Explain the difference between present_calibrations and matched_calibrations
#' @details The function takes pairs of taxon names in a secondary calibrations data frame,
#' and looks for them in the vector of tip labels of the tree. If both are present,
#' then it gets the node that represents the Most Recent
#' Common Ancestor for that pair of taxa in the tree.
#' @export
match_all_calibrations <- function(phy, calibrations) {
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
  calibrations <- calibrations[which(calibrations$taxonA %in% phy$tip.label), ]
  calibrations <- calibrations[which(calibrations$taxonB %in% phy$tip.label), ]
  if (nrow(calibrations) == 0) {
    warning("Taxon name pairs in 'calibrations' do not match 'phy' tip labels.")
    return(NA)
  }
  age_column <- grep("age", tolower(names(calibrations))) # get columns with "age" in their names
  target_tree_nodes <- sapply(
    seq(nrow(calibrations)),
    function(i) {
      phytools::findMRCA(
        tree = phy, type = "node",
        tips = as.character(calibrations[i, c("taxonA", "taxonB")])
      )
    }
  )
  target_tree_nodes <- target_tree_nodes - ape::Ntip(phy)
  all_nodes <- sort(unique(target_tree_nodes))
  # get the node age distribution (ages taken from the first column with "age" in its name):
  all_ages <- lapply(all_nodes, function(i) calibrations[target_tree_nodes == i, age_column[1]])
  # any(sapply(all_ages, is.null)) # if FALSE, all nodes have at least one calibration.
  rowsies <- !duplicated(target_tree_nodes)
  target_tree_nodes2 <- target_tree_nodes[rowsies]
  calibrations2 <- calibrations[rowsies, ]
  calibrations2 <- calibrations2[order(target_tree_nodes2), ]
  calibrations2$MaxAge <- sapply(all_ages, max)
  calibrations2$MinAge <- sapply(all_ages, min)
  calibrations2$NodeNames <- paste0("n", all_nodes)
  if (all(all_nodes < ape::Ntip(phy))) {
    all_nodes_numbers <- all_nodes + ape::Ntip(phy)
    node_index <- "consecutive"
  } else {
    all_nodes_numbers <- all_nodes
    node_index <- "node_number"
  }
  # so we can rename all nodes of interest to match our labels, make sure node.label is null:
  phy$node.label <- NULL
  # all nodes need to be named so that make_bladj_tree runs properly:
  phy <- tree_add_nodelabels(
    tree = phy,
    node_index = node_index
  )
  phy$calibration_distribution <- stats::setNames(
    all_ages,
    all_nodes_numbers
  )
  return(list(
    phy = phy,
    matched_calibrations = calibrations2,
    present_calibrations = calibrations
  ))
}
