#' Use the BLADJ algorithm to get a chronogram from a tree topology for which you have age data for some of its nodes.
#'
#' @description The function takes a tree topology and uses the BLADJ algorithm
#' implemented with [phylocomr::ph_bladj()] to assign node ages and branch lengths, given a
#' set of fixed node ages and respective node names.
#' @param nodenames A character vector with names of nodes in tree with known ages
#' @param nodeages A numeric vector with the actual ages of named nodes
#' @inheritParams tree_fix_brlen
#' @return A `phylo` object.
#' @details Input `tree` can be dated or not, `$edge.length` is ignored.
#' Ages given in `nodeages` are fixed on their corresponding nodes given in `nodenames`.
#' @export
make_bladj_tree <- function(tree = NULL,
                            nodenames = NULL,
                            nodeages = NULL) {
  ############################################################################
  # initial checks
  ############################################################################
  phy <- tree_check(tree = tree, dated = FALSE)
  # needs to be fully resolved to work with bladj?
  # ape::is.binary(phy)
  if (is.null(phy$node.label)) {
    stop("'tree' must have node labels")
  }
  if (!is.null(phy$edge.length)) {
    phy$edge.length <- NULL
  }
  if (!is.character(nodenames)) {
    stop("'nodenames' must be a character vector")
  }
  if (!is.numeric(nodeages)) {
    stop("'nodeages' must be a numeric vector")
  }
  # phylocomr takes nodenames with spaces in it as two items, so it fails to read
  # the table of ages later, thinking that it has 3 columns instead of 2
  nodenames <- gsub(" ", "_", nodenames)
  phy$node.label <- gsub(" ", "_", phy$node.label)
  m <- match(nodenames, phy$node.label)
  if (any(is.na(m))) {
    # say which nodenames are not in phy$node.label
    warning("Not all 'nodenames' are in 'phy$node.label'; these will be ignored:\n",
            nodenames[is.na(m)])
    # remove the node ages and node names that are not in phy$node.label:
    nodenames <- nodenames[!is.na(m)]
    nodeages <- nodeages[!is.na(m)]
  }
  if (length(nodenames) != length(nodeages)) {
    stop("'nodenames' and 'nodeages' must have the same length")
  }
  ############################################################################
  # create the dataframe that goes into bladj:
  ############################################################################
  ages_df <- data.frame(a = nodenames,
                        b = nodeages)
  # to be used by ph_bladj, phy cannot have more classes
  class(phy) <- "phylo"
  ############################################################################
  # run bladj
  ############################################################################
  new.phy <- phylocomr::ph_bladj(ages = ages_df, phylo = phy)
  attributes(new.phy) <- NULL
  # new.phy <- phytools::read.newick(text = new.phy) # usually more robust
  new.phy <- ape::read.tree(text = new.phy)
  # to keep the same names as original phy (bladj modifies all names to lowercase):
  ############################################################################
  # format the chronogram
  ############################################################################
  phy <- phylo_tiplabel_space_to_underscore(phy)
  phy$tip.label <- gsub(":", "", phy$tip.label) # one tip label in Hedges et al. 2015 chronogram
  phy$tip.label <- gsub("\\(", "-", phy$tip.label) # one tip label in Jetz et al. 2012 chronogram # "bernieria_madagascariensis_(gmelin,_1789)"
  phy$tip.label <- gsub(")", "", phy$tip.label) # one tip label in Jetz et al. 2012 chronogram # "bernieria_madagascariensis_(gmelin,_1789)"
  phy$tip.label <- gsub(",", "", phy$tip.label) # one tip label in Jetz et al. 2012 chronogram # "bernieria_madagascariensis_(gmelin,_1789)"
  # new.phy$tip.label <- gsub("-", "", new.phy$tip.label) # one tip label in Jetz et al. 2012 chronogram # "bernieria_madagascariensis_(gmelin,_1789)"
  index <- match(tolower(phy$tip.label), new.phy$tip.label)
  # index2 <- match(new.phy$tip.label, tolower(phy$tip.label))
  # tolower(new.phy$tip.label)[is.na(index2)]
  # tolower(phy$tip.label)[is.na(index)]
  # new.phy$tip.label[is.na(index2)]
  new.phy$tip.label[index] <- phy$tip.label
  # Error in new.phy$tip.label[index] <- phy$tip.label :
  # NAs are not allowed in subscripted assignments
  return(new.phy)
}
# new.phy_save <- new.phy
