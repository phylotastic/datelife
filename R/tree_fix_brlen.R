# functions to fix negative or null branch lengths from an ultranetric treem, and to add dates of missing taxa at random

#' Take a tree with branch lengths and fix negative or zero length branches.
#' @param tree A tree either as a newick character string or as a `phylo` object.
#' @param fixing_criterion A character vector specifying the type of branch length to be fixed: "negative" or "zero" (the number 0 is also allowed).
#' @param fixing_method A character vector specifying the method to fix branch lengths: "bladj", "mrbayes" or a number to be assigned to all branches meeting fixing_criterion
#' @param ultrametric Boolean indicating whether to force ultrametric or not.
#' @return A `phylo` object with no negative or zero branch lengths.
#' @export
tree_fix_brlen <- function(tree = NULL, fixing_criterion = "negative", fixing_method = 0, ultrametric = TRUE) {
  phy <- tree_check(tree = tree, brlen = TRUE, dated = FALSE)
  # phytools::force.ultrametric fixes negative branch lengths at some extent, but it sometimes leaves tiny neg branches still; those can be fixed later in here.
  # but it is helpful to keep the tree ultrametric, that is implememted within force_ultrametric:
  if (ultrametric) {
    phy <- force_ultrametric(phy)
  }
  # enhance: allow to fix both neg and zero at the same time!
  fixing_criterion <- match.arg(
    arg = as.character(fixing_criterion), choices =
      c("negative", "zero", "0"), several.ok = FALSE
  )
  if (fixing_criterion == "negative") {
    index <- which(phy$edge.length < 0) # identifies edge numbers with negative edge lengths value
  } else {
    index <- which(phy$edge.length == 0) # identifies edge numbers with null/zero edge lengths value
  }
  if (length(index) == 0) {
    return(phy)
  }
  if (!is.numeric(fixing_method)) {
    fixing_method <- match.arg(fixing_method, c("bladj", "mrbayes"))
  } else { # chunk for neg or zero br len to zero or any number determined by user
    # stop()
    for (i in index) {
      # snode <- pos.phy$edge[i,1]
      # pool  <- pos.phy$edge[seq(nrow(pos.phy$edge))[-i], 1]
      # sisedge <- which(pool==snode) # determines position of sister edge
      # pos.phy$edge.length[sisedge] <- pos.phy$edge.length[sisedge] - pos.phy$edge.length[i]
      # adds neg branch length to sister branch, should add error to both sides???? or only to the daughter branches??
      cnode <- phy$edge[i, 2]
      dauedge <- which(phy$edge[, 1] == cnode)
      phy$edge.length[dauedge] <- phy$edge.length[dauedge] + phy$edge.length[i] + fixing_method[1]
      phy$edge.length[i] <- fixing_method[1]
      fixed.phy <- phy
    }
  }

  if (any(fixing_method == c("bladj", "mrbayes"))) { # chunk for bladj and mrbayes
    phy <- tree_add_nodelabels(tree = phy) # all nodes need to be named
    cnode <- phy$edge[index, 2] # we assume that the negative edge length is the one that needs to be changed (but it could be the sister edge that should be shorter)
    tofix <- cnode - length(phy$tip.label) # so, we use the node number of the crown node of the negative branch lengths
    if (fixing_method == "bladj") {
      fixed.phy <- make_bladj_tree(
        tree = phy,
        nodenames = phy$node.label[-tofix],
        nodeages = tree_get_node_data(
          tree = phy,
          node_data = "node_age"
        )$node_age[-tofix]
      )
    }
    if (fixing_method == "mrbayes") {
      mrbayes.file <- paste0(deparse(substitute(tree)), "_brlen_fixed.nexus") # make "phylo" an argument?
      ncalibration <- tree_get_node_data(
        tree = phy,
        node_data = c("descendant_tips_label", "node_age")
      )
      ncalibration <- lapply(ncalibration, "[", seq(phy$Nnode)[-tofix])
      phy <- tree_add_outgroup(tree = phy, outgroup = "fake_outgroup")
      fixed.phy <- make_mrbayes_tree(
        constraint = phy, taxa = phy$tip.label,
        ncalibration = ncalibration, age_distribution = "fixed",
        root_calibration = FALSE,
        missing_taxa = NULL, mrbayes_output_file = mrbayes.file
      )
      fixed.phy <- ape::drop.tip(fixed.phy, "fake_outgroup")
    }
  }
  return(fixed.phy)
}
