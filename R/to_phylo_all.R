#' Gets all ages per taxon pair from a distance matrix
#' Internal function used in summary_matrix_to_phylo_all().
#' @return A `data.frame` of pairwise ages, with row number equal to the combinatory
#' of column names (or row names), estimated as `ncol(summ_matrix)^2 - sum(1:(ncol(summ_matrix)-1))`.
#' @inheritParams summary_matrix_to_phylo
summarize_summary_matrix <- function(summ_matrix) {
  ############################################################################
  ############################################################################
  # to compute the final length of the data frame do: ncol(xx)^2 - sum(1:(ncol(xx)-1))
  # calibrations <- matrix(nrow = ncol(xx)^2 - sum(1:(ncol(xx)-1)), ncol = 3)
  # identify if SDM matrix has some negative values; extract taxon names:
  ############################################################################
  ############################################################################
  negs <- which(summ_matrix < 0)
  neg_names <- rownames(summ_matrix)[ceiling(negs / nrow(summ_matrix))]
  ############################################################################
  ############################################################################
  # extract pairwisee ages from summ_matrix:
  ############################################################################
  ############################################################################
  ages <- tA <- tB <- c()
  for (i in seq(ncol(summ_matrix))) {
    ages <- c(ages, summ_matrix[1:i, i])
    tA <- c(tA, rownames(summ_matrix)[1:i])
    tB <- c(tB, rep(colnames(summ_matrix)[i], i))
  }
  tA <- gsub(" ", "_", tA)
  tB <- gsub(" ", "_", tB)
  ############################################################################
  ############################################################################
  # store pairwise ages as data frame:
  ############################################################################
  ############################################################################
  calibrations <- data.frame(Age = ages, taxonA = tA, taxonB = tB, stringsAsFactors = FALSE)
  calibrations <- calibrations[!is.na(calibrations[, "Age"]), ] # get rid of NaN and NAs
  calibrations <- calibrations[calibrations[, "Age"] != 0, ] # get rid of 0's
  calibrations <- calibrations[calibrations[, "Age"] > 0, ] # get rid of negative values too
  if (any(is.na(calibrations[, "Age"]))) {
    warning("for some reason there are still NAs in the matrix")
  }
  # SDM summary matrix sometimes has negative values,
  # bc ages are transformed to be approximated in a similar way as a linear regression
  ############################################################################
  ############################################################################
  # Order calibrations by taxon name and age, and return
  ############################################################################
  ############################################################################
  calibrations <- calibrations[order(calibrations$taxonA,
                                     calibrations$Age,
                                     calibrations$taxonB),]

  rownames(calibrations) <- NULL
  return(calibrations)
}

#' Get minimum, median, mean, midpoint, and maximum summary chronograms from a
#' summary matrix of a `datelifeResult` object.
#' @inheritParams summary_matrix_to_phylo
#' @inheritDotParams get_otol_synthetic_tree
#' @return A `multiPhylo` object of length 5. It contains min, mean, median, midpoint, and max summary chronograms.
#' @details
#' With this function users can choose the minimum, mean or maximum ages from
#' the summary matrix as calibration points to get a single summary chronogram.
#' Users get all three summary chronograms in a `multiPhylo` object.
# Modified from `get_all_summaries()` function in `data-raw/datelife_examples.R`
#' @export
summary_matrix_to_phylo_all <- function(summ_matrix,
                                        datelife_query = NULL,
                                        target_tree = NULL,
                                        total_distance = TRUE,
                                        ...) {

  if (!inherits(summ_matrix, "matrix") & !inherits(summ_matrix, "data.frame")) {
    message("'summ_matrix' argument is not a matrix")
    return(NA)
  }
  if (!is.null(datelife_query)) {
    input_ott_match <- suppressMessages(check_ott_input(input = datelife_query))
    # match inputt_ott_match and unique(c(colnames(summ_matrix), rownames(summ_matrix)))
    # change the names in target tree to the names from summ_matrix (which are the ones that come from the original query)
  }
  # summ_matrix <- data.frame(summ_matrix)
  # everything up to patristic_matrix_to_phylo ok if it is a data frame too
  if (inherits(summ_matrix, "data.frame")) {
    summ_matrix <- as.matrix(summ_matrix)
    colnames(summ_matrix) <- gsub("\\.", " ", colnames(summ_matrix))
  }
  if (total_distance) {
    summ_matrix <- summ_matrix * 0.5 # bc it's total distance tip to tip
  }
  ############################################################################
  ############################################################################
  # get a backbone tree topology if none is provided
  ############################################################################
  ############################################################################
  if (!inherits(target_tree, "phylo")) {
    message("... No target_tree was provided, obtaining a tree topology from Open Tree synthetic phylogeny.")
    target_tree <- suppressMessages(get_otol_synthetic_tree(input = colnames(summ_matrix), ...))
    if (!inherits(target_tree, "phylo")) {
      # enhance: we should find a better way to do this, but it should be ok for now:
      message(paste("Obtaining a topology from Open Tree failed."))
      message(paste("... Converting patristic matrix to phylo."))
      target_tree <- suppressWarnings(suppressMessages(patristic_matrix_to_phylo(summ_matrix, ultrametric = TRUE)))
      # target_tree <- consensus(phyloall, p = 0.5) # can't use consensus here: bc not all trees have the same number of tips
    }
    target_tree <- ape::collapse.singles(target_tree)
  }
  if (!inherits(target_tree, "phylo")) {
    message("'target_tree' is missing or not a 'phylo' object and a backbone tree could not be constructed; returning 'NA'")
    message("Hint: Was summ_matrix constructed from an object with no good groves? Try running 'get_best_grove' first.")
    # enhance: add a more formal test of best grove
    return(NA)
  }
  target_tree$edge.length <- NULL
  target_tree$edge.length.original <- NULL
  target_tree$tip.label <- gsub(" ", "_", target_tree$tip.label)
  ############################################################################
  ############################################################################
  # test that taxa in matrix are all in target tree tip labels
  ############################################################################
  ############################################################################
  rownames(summ_matrix) <- gsub(" ", "_", rownames(summ_matrix))
  colnames(summ_matrix) <- gsub(" ", "_", colnames(summ_matrix))
  # find taxa missing in target tree and remove them from summ_matrix
  missing <- is.na(match(colnames(summ_matrix), target_tree$tip.label))
  whichmiss <- colnames(summ_matrix)[missing]
  if (any(missing)) {
    message("Some taxa in summ_matrix are not in target_tree (", paste0(whichmiss, collapse = ", "), ")")
    missingrow <- is.na(match(rownames(summ_matrix), target_tree$tip.label))
    summ_matrix <- summ_matrix[!missingrow, !missing]
  }
  # to be get_all_calibrations.data.frame:
  ############################################################################
  ############################################################################
  # extract pairwise ages
  ############################################################################
  ############################################################################
  calibrations <- summarize_summary_matrix(summ_matrix)
  ##############################################################################
  ##############################################################################
  # get mrca node numbers from calibration pairs
  ############################################################################
  ############################################################################
  # start of use_all_calibrations_bladj, that contains match_all_calibrations
  # something like:
  # use_all_calibrations_bladj(phy = target_tree, calibrations = calibrations, type = use)
  # start of match_all_calibrations or mrca calibrations:
  # get the coincident node numbers:
  # ape::is.binary(target_tree)
  target_tree_nodes <- sapply(seq(nrow(calibrations)), function(i) {
    phytools::findMRCA(
      tree = target_tree,
      tips = as.character(calibrations[i, c("taxonA", "taxonB")]),
      type = "node"
    )
  })
  # Start node numbers at 1:
  target_tree_nodes <- target_tree_nodes - ape::Ntip(target_tree)
  ##############################################################################
  ##############################################################################
  # get ages per target node (node age distribution)
  ############################################################################
  ############################################################################
  all_nodes <- sort(unique(target_tree_nodes))
  all_ages <- lapply(all_nodes, function(i) calibrations[target_tree_nodes == i, "Age"])
  # any(sapply(all_ages, is.null)) # if FALSE, all nodes have at least one calibration.
  calibrations2 <- data.frame(MRCA = paste0("n", all_nodes),
                              MinAge = sapply(all_ages, min),
                              MaxAge = sapply(all_ages, max))
  # calibrations2$MRCA is a factor so have to be made as.character to work with bladj
  ##############################################################################
  ##############################################################################
  # rename nodes on target tree and age distributions:
  ############################################################################
  ############################################################################
  # give node names to the distribution of node ages
  names(all_ages) <- paste0("n", all_nodes)
  # then make sure node labels are null,
  # so we can rename all nodes of interest to match our labels
  target_tree$node.label <- NULL
  # all nodes need to be named so make_bladj_tree runs properly:
  node_index <- ifelse(all(all_nodes < ape::Ntip(target_tree)),
                       "from_1",
                       "node_number")
  # if (all(all_nodes < ape::Ntip(target_tree))) {
  #   # all_nodes_numbers <- all_nodes + ape::Ntip(target_tree)
  #   node_index <- "from_1"
  # } else {
  #   # all_nodes_numbers <- all_nodes
  #   node_index <- "node_number"
  # }
  target_tree <- tree_add_nodelabels(tree = target_tree, node_index = node_index)
  # end of match_all_calibrations or mrca calibrations
  ##############################################################################
  ##############################################################################
  # use bladj to date the tree using node age distributions:
  ############################################################################
  ############################################################################
  node_ages_midpoint <- sapply(seq(nrow(calibrations2)),
                        function(i) sum(calibrations2[i, c("MinAge", "MaxAge")]) / 2)
  new_phy_midpoint <- make_bladj_tree(
    tree = target_tree, nodenames = as.character(names(all_ages)),
    nodeages = node_ages_midpoint
  )
  new_phy_min <- make_bladj_tree(
    tree = target_tree, nodenames = as.character(names(all_ages)),
    nodeages = sapply(all_ages, min)
  )
  new_phy_max <- make_bladj_tree(
    tree = target_tree, nodenames = as.character(names(all_ages)),
    nodeages = sapply(all_ages, max)
  )
  new_phy_mean <- make_bladj_tree(
    tree = target_tree, nodenames = as.character(names(all_ages)),
    nodeages = sapply(all_ages, mean)
  )
  new_phy_median <- make_bladj_tree(
    tree = target_tree, nodenames = as.character(names(all_ages)),
    nodeages = sapply(all_ages, stats::median)
  )
  # end use_all_calibrations_bladj
  res <- c(new_phy_min, new_phy_max, new_phy_median, new_phy_mean, new_phy_midpoint)
  attributes(res)$node_age_distributions <- all_ages
  names(res) <- c("min", "max", "median", "mean", "midpoint")
  class(res) <- "multiPhylo"

  return(res)
}
