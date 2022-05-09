#' Match calibrations to nodes of a given tree
#'
#' @description \code{match_all_calibrations} searches a given tree for the most recent common
#'   ancestor (mrca) of all taxon name pairs in a `datelifeCalibration`. It uses [phytools::findMRCA()].
#'
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param calibrations A `calibrations` object, an output of
#'   [extract_calibrations_phylo()].
#' @return A list of two elements:
#' \describe{
#' 	\item{phy}{A `phylo` object with nodes renamed with [tree_add_nodelabels()].}
#' 	\item{matched_calibrations}{A `matchedCalibrations` object, which is the input `calibrations`
#'    object with two additional columns storing results from the mrca search with
#'    [phytools::findMRCA()]: `$mrca_node_number` and `$mrca_node_name`.}
#' 	}
#' @details The function takes pairs of taxon names in a secondary calibrations data frame,
#' and looks for them in the vector of tip labels of the tree. If both are present,
#' then it gets the node that represents the most recent
#' common ancestor (mrca) for that pair of taxa in the tree.
#' Nodes of input `phy` can be named or not.
# #' @export
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

# attempt to simplify the code to get mrca_nodes:
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
  return(list(phy = phy,
              matched_calibrations = structure(calibrations,
                                               class = c("data.frame",
                                                         "matchedCalibrations")
                                     )
        )
  )
}
#' Summarize a `matchedCalibrations` object
#' `summary.matchedCalibrations` gets the node age distribution from a `matchedCalibrations` object.
#' @param object A `matchedCalibrations` object, usually an element of the output of [match_all_calibrations()].
#' @param ... Further arguments passed to or from other methods.
#' @return A `summaryMatchedCalibrations` object, which is a list of two `matchedCalibrations` objects:
#' \describe{
#' 	\item{not_in_phy}{A `data.frame` subset of input `matchedCalibrations` object
#' 	  containing taxon name pairs that were not present in the given tree. `NULL`
#' 	  if all input taxon names are found in the given tree.}
#' 	\item{in_phy}{A `data.frame` subset of input `matchedCalibrations` object
#' 	  containing all taxon name pairs that were present in the given tree.}
#' }
#' @details Columns `in_phy$mrca_node_name` and `in_phy$reference` are factors.
#' @export
summary.matchedCalibrations <- function(object, ...) {
  all_nodes <- sort(unique(object$mrca_node_number))

  # Subset the data.frame:
  not_in_phy_rows <- which(is.na(object$mrca_node_number))
  message1 <- c()
  if (length(not_in_phy_rows) > 0) {
    not_in_phy <- object[not_in_phy_rows, ]
    in_phy <- object[-not_in_phy_rows, ]
    message1 <- c(message1, "Not all taxon name pairs are in 'phy'.")
  } else {
    message1 <- c(message1, "All taxon name pairs are in 'phy'.")
    not_in_phy <- NULL
    in_phy <- object
  }
  in_phy$mrca_node_name <- as.factor(in_phy$mrca_node_name)
  in_phy$reference <- as.factor(in_phy$reference)
  # is MaxAge and MinAge the same value?
  if (all(in_phy$MaxAge == in_phy$MinAge)) {
    message1 <- c(message1,
                  "\n'MaxAge' and 'MinAge' columns in input 'matchedCalibrations' /
                   have the same values.")
  }
  message("Success!")
  message(message1)
  return(structure(list(not_in_phy = not_in_phy, in_phy = in_phy),
                   class = c("list", "summaryMatchedCalibrations")))
}


# # Get the node age distribution:
# # So we take data from the first column name matching the word "age":
# age_column <- grep("age", tolower(names(object))) # get columns with "age" in their names
# # create a list with node age distributions:
# ages_distribution <- lapply(all_nodes, function(i) object[all_nodes == i, age_column[1]])
#
# # get the references for ages in ages_distribution:
# ages_references <- lapply(all_nodes, function(i) calibrations[all_nodes == i, "reference"])
#
#
# # any(sapply(ages_distribution, is.null)) # if FALSE, all nodes have at least one calibration.
# rowsies <- !duplicated(mrca_nodes)
# mrca_nodes2 <- mrca_nodes[rowsies]
# calibrations2 <- calibrations[rowsies, ]
# calibrations2 <- calibrations2[order(mrca_nodes2), ]
# calibrations2$MaxAge <- sapply(ages_distribution, max)
# calibrations2$MinAge <- sapply(ages_distribution, min)
# calibrations2$NodeNames <- paste0("n", all_nodes)
