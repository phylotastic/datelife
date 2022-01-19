# functions to process a tree in newick or phylo format, and get information on its properties
# most of them are useful to process a tree so it can be taken up by bladj or mrbayes

#' Get node numbers, node names, descendant tip numbers and labels of nodes from any tree, and node ages from dated trees.
#' @inheritParams tree_fix_brlen
#' @param node_data A character vector containing one or all from: "node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label"
#' @param nodes Numeric vector with node numbers from which you want to obtain data. Default to NULL: obtain data for all nodes in the tree.
#' @return A list
#' @export
tree_get_node_data <- function(tree = NULL, nodes = NULL, node_data = c("node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label")) {
  node_data <- match.arg(arg = node_data, choices = c("node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label"), several.ok = TRUE)
  if ("node_age" %in% node_data) {
    phy <- tree_check(tree = tree, dated = TRUE)
  } else {
    phy <- tree_check(tree = tree)
  }
  if (!is.numeric(nodes)) {
    nn <- phylo_get_node_numbers(phy)
  } else {
    nn <- nodes
  }
  res <- vector(mode = "list")
  if ("node_number" %in% node_data) {
    res <- c(res, list(node_number = nn))
  }
  if ("node_label" %in% node_data) {
    res <- c(res, list(node_label = phy$node.label))
  }
  if ("node_age" %in% node_data) {
    res <- c(res, list(node_age = ape::branching.times(phy)))
  }
  if (any(c("descendant_tips_number", "descendant_tips_label") %in% node_data)) {
    dt_num <- lapply(nn, function(x) tree_node_tips(tree = phy, node = x)) # tip numbers stemming from node
    names(dt_num) <- nn
    dt_lab <- lapply(dt_num, function(x) phy$tip.label[x]) # tip labels corresponding to those tips numbers
    names(dt_lab) <- nn
    if ("descendant_tips_number" %in% node_data) {
      res <- c(res, list(descendant_tips_number = dt_num))
    }
    if ("descendant_tips_label" %in% node_data) {
      res <- c(res, list(descendant_tips_label = dt_lab))
    }
  }
  return(res)
}

#' Adds labels to nodes with no assigned label
#' @inheritParams tree_fix_brlen
#' @param node_prefix Character vector. If length 1, it will be used to name all nodes with no labels, followed by a number which can be the node_number or consecutive, as specified in node_index.
#' @param node_index Character vector. Choose between "consecutive" and "node_number" as numeric index for node labels. It will use consecutive numbers from 1 to total node number in the first case and phylo node numbers in the second case (i.e, from Ntip + 1).
#' @return A phylo object
#' @export
tree_add_nodelabels <- function(tree = NULL, node_prefix = "n", node_index = "node_number") {
  phy <- tree_check(tree = tree)
  node_index <- match.arg(arg = node_index, choices = c("consecutive", "node_number"), several.ok = FALSE)
  if ("node_number" %in% node_index) {
    node_number <- phylo_get_node_numbers(phy = phy)
  }
  if ("consecutive" %in% node_index) {
    node_number <- seq(phy$Nnode)
  }
  if (is.null(phy$node.label)) {
    phy$node.label <- paste0(node_prefix, node_number)
  } else {
    en <- which(phy$node.label == "")
    phy$node.label[en] <- paste0(node_prefix, node_number[en])
  }
  return(phy)
}

#' Gets node numbers from any phylogeny
#' @inheritParams phylo_check
#' @return A numeric vector with node numbers
#' @export
phylo_get_node_numbers <- function(phy) {
  # node_numbers <- (length(phy$tip.label) + 1):(length(phy$tip.label) + phy$Nnode)
  node_numbers <- unique(phy$edge[, 1][which(phy$edge[, 1] > length(phy$tip.label))])
  node_numbers <- node_numbers[order(node_numbers)]
  return(node_numbers)
}

#' Function to add an outgroup to any phylogeny, in phylo or newick format
#' @inheritParams tree_fix_brlen
#' @param outgroup A character vector with the name of the outgroup. If it has length>1, only first element will be used.
#' @return A phylo object with no root edge.
#' @export
tree_add_outgroup <- function(tree = NULL, outgroup = "outgroup") {
  phy <- tree_check(tree = tree)
  outgroup_edge <- c()
  ingroup_edge <- c()
  root_edge <- c() # in case we want to allow users to add a root edge
  if (!is.null(phy$edge.length)) {
    mbt <- max(ape::branching.times(phy))
    if (!is.null(phy$root.edge)) {
      # root_edge <- paste0(":", phy$root.edge)
      rr <- phy$root.edge
      phy <- phy[-which(names(phy) == "root.edge")]
      class(phy) <- "phylo"
    } else {
      rr <- mbt * 0.10
    }
    ingroup_edge <- paste0(":", rr)
    outgroup_edge <- paste0(":", mbt + rr)
  }
  phy <- gsub(";", "", ape::write.tree(phy))
  phy <- paste("(", phy, ingroup_edge, ",", outgroup[1], outgroup_edge, ")", root_edge, ";", sep = "")
  # phy <-  phytools::read.newick(text = phy)
  phy <- ape::read.tree(text = phy) # tree looses its root length anyways when read by read.tree or read.newick

  return(phy)
}

#' Checks if `phy` is a `phylo` object and/or a chronogram.
#' @param phy A `phylo` object.
#' @param brlen Boolean. If `TRUE` it checks if `phylo` object has branch lengths.
#' @param dated Boolean. If `TRUE` it checks if `phylo` object is ultrametric.
#' @return Nothing
#' @export
phylo_check <- function(phy = NULL, brlen = FALSE, dated = FALSE) {
  if (!inherits(phy, "phylo")) {
    stop("tree is not a phylo object")
  }
  if (brlen) {
    if (!phylo_has_brlen(phy = phy)) {
      stop("tree must have branch lengths")
    }
  }
  if (dated) {
    if (!ape::is.ultrametric(phy, option = 2)) {
      warning("branch lengths in tree should be proportional to time")
      stop("tree must be ultrametric") # really?  # Think how to incorporate trees with extinct taxa
    }
  }
}

#' Checks if a tree is a phylo class object otherwise it uses input_process. Additionally it can check if tree is a chronogram with phylo_check
#' @inheritParams tree_fix_brlen
#' @inheritDotParams phylo_check -phy
#' @return If tree is correctly formatted, it returns a `phylo` object.
#' @export
tree_check <- function(tree = NULL, ...) {
  if (!inherits(tree, "phylo")) {
    tree <- input_process(input = tree)
  }
  phylo_check(phy = tree, ...)
  return(tree)
}

#' To get tip numbers descending from any given node of a tree
#' @inheritParams phytools::getDescendants
#' @inheritParams tree_fix_brlen
#' @return A numeric vector with tip numbers descending from a node
#' @export
tree_node_tips <- function(tree = NULL, node = NULL, curr = NULL) {
  phy <- tree_check(tree = tree, dated = FALSE)
  des <- phytools::getDescendants(tree = phy, node = node, curr = NULL)
  tips <- des[which(des <= length(phy$tip.label))]
  # we could just use ape::prop.part here -replacing last two lines, it might be faster and using less memory
  # ape::prop.part gets tip numbers in node order
  # we would just need to exclude data from nodes we don't want
  return(tips)
}
