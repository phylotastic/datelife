#' Add missing taxa to a dated tree and fabricate node ages for these missing taxa.
#'
#' @description This function adds missing taxa to a chronogram given in `dated_tree`.
#' It is still work in progress.
#' @param dated_tree a tree (newick or phylo) with branch lengths proportional to absolute time
#' @inheritParams  missing_taxa_check
#' @param dating_method The method used for tree dating, options are "mrbayes" and "bladj".
#' @param adding_criterion Only used when `dating_method = "mrbayes"`. A character vector to specify how missing_taxa should be added to dated_tree.
#' 	 Choose one of:
#' \describe{
#' 	\item{adding_method = "random"}{missing_taxa will be added at random to dated_tree.
#' 	}
#' 	\item{adding_method = "taxonomy"}{taxa will be added to dated_tree following a dataframe with taxonomic assignations given in missing_taxa argument. If no dataframe is given, OpenTree's reference taxonomy will be used.
#' 	}
#' 	\item{adding_method = "tree"}{taxa will be added to dated_tree following a tree given in missing_taxa argument. If no tree is given, OpenTree's synthetic tree will be used.
#' 	}
#' }
#' @inheritParams make_mrbayes_runfile
#' @return A `phylo` object.
tree_add_dates <- function(dated_tree = NULL,
                           missing_taxa = NULL,
                           dating_method = "mrbayes",
                           adding_criterion = "random",
                           mrbayes_output_file = "mrbayes_tree_add_dates.nexus") {
  dated_tree <- tree_check(tree = dated_tree, dated = TRUE)
  missing_taxa <- missing_taxa_check(missing_taxa = missing_taxa, dated_tree = dated_tree)
  dating_method <- match.arg(dating_method, c("bladj", "mrbayes"))
  adding_criterion <- tryCatch(match.arg(adding_criterion, c("random", "taxonomy", "tree")), error = function(e) "random") # if it does not match any it is assigned to NULL
  if (dating_method == "bladj") {
    # we need to add a missing_taxa check here. We can only use bladj if missing_taxa is a tree
    # and fully resolved
    if (inherits(missing_taxa, "phylo")) {
      missing_taxa_phy <- missing_taxa
    } else {
      if (is.data.frame(missing_taxa)) {
        all_taxa <- unique(c(dated_tree$tip.label, levels(missing_taxa$taxon)))
      }
      if (is.vector(missing_taxa)) {
        all_taxa <- unique(c(dated_tree$tip.label, missing_taxa))
      }
      missing_taxa_phy <- get_otol_synthetic_tree(input = all_taxa) # tip labes have underscores already
      # this does not always recovers all taxa missing
      # add a warning and a suggestion to rerun tree_add_dates with remaining absent taxa
      # add absent_taxa element here too...
    }
    # enhance: check the following warning in congruify.phylo
    # In if (class(stock) == "phylo") { :
    #   the condition has length > 1 and only the first element will be used)
    constraint_tree <- suppressWarnings(geiger::congruify.phylo(
      reference = phylo_tiplabel_space_to_underscore(dated_tree),
      target = phylo_tiplabel_space_to_underscore(missing_taxa_phy), scale = NA,
      ncores = 1
    ))
    # add congruified nodes as node labels to missing_taxa_phy
    dated_tree_nodes <- sapply(seq(nrow(constraint_tree$calibrations)), function(i) {
      phytools::findMRCA(
        tree = constraint_tree$target,
        tips = as.character(constraint_tree$calibrations[i, c("taxonA", "taxonB")]),
        type = "node"
      )
    })
    dated_tree_nodes <- dated_tree_nodes - ape::Ntip(constraint_tree$target)
    # missing_taxa_phy$node.label[dated_tree_nodes] <- constraint_tree$calibrations$MRCA
    # cannot use hash number to name nodes, bladj collapses. So using "congNumber"
    missing_taxa_phy$node.label[dated_tree_nodes] <- paste0("cong", seq(nrow(constraint_tree$calibrations)))
    missing_taxa_phy <- tree_add_nodelabels(tree = missing_taxa_phy) # all nodes need to be named so make_bladj_tree runs properly
    # this adds random names to unnamed nodes, but they have to coincide between target and reference
    # new.phy <- make_bladj_tree(tree = missing_taxa, nodenames = dated_tree$node.label, nodeages = tree_get_node_data(tree = dated_tree, node_data = "node_age")$node_age)
    new.phy <- make_bladj_tree(tree = missing_taxa_phy, nodenames = missing_taxa_phy$node.label[dated_tree_nodes], nodeages = sapply(seq(nrow(constraint_tree$calibrations)), function(i) sum(constraint_tree$calibrations[i, c("MinAge", "MaxAge")]) / 2))
  }
  if (dating_method == "mrbayes") {
    dated_tree <- tree_add_outgroup(tree = dated_tree, outgroup = "an_outgroup") # we need to add a fake outgroup, otherwise mrbayes won't respect the root age
    ncalibration <- tree_get_node_data(tree = dated_tree, node_data = c("node_age", "descendant_tips_label"))
    # we need to be more specific in the way it uses missing taxa next. If it is a tree, then missing_taxa goes as the constraint. If it is a vector, then it just goes as missing taxa. If it is a data frame, we should call pastis.
    new.phy <- make_mrbayes_tree(constraint = dated_tree, ncalibration = ncalibration, missing_taxa = missing_taxa, mrbayes_output_file = mrbayes_output_file)
    new.phy <- ape::drop.tip(new.phy, "an_outgroup")
  }
  return(new.phy)
}


#' Checks that missing_taxa argument is ok to be used by make_mrbayes_runfile inside tree_add_dates functions.
#' @param  missing_taxa A tree, a data frame or a vector enlisting all missing taxa you want to include.
#' \describe{
#' 	\item{A tree}{Either as a phylo object or as a newick character string.
#' 		It contains all taxa that you want at the end, both missing and non missing.
#' 		This tree will be used as a hard constraint.}
#' 	\item{A `data.frame`}{It contains two columns named "taxon" and "clade".
#' 		The first one contains a character vector of missing taxon names.
#' 		The second one contains a character or numeric vector of nodes from a
#'    constraint tree to which each taxon will be assigned.}
#' 	\item{A character vector}{It contains the names of the missing taxa.
#' 		They will be added at random to the constraint tree.}
#' }
#' @inheritParams tree_add_dates
#' @return A phylo object, a newick character string or a dataframe with taxonomic assignations
#' @export
missing_taxa_check <- function(missing_taxa = NULL, dated_tree = NULL) {
  badformat <- TRUE
  if (is.data.frame(missing_taxa)) { # or is.matrix??
    if ("taxon" %in% names(missing_taxa)) {
      badformat <- FALSE
      if (length(missing_taxa) > 1 & !"clade" %in% names(missing_taxa)) {
        badformat <- TRUE
      }
    }
    # check PastisData format to further check data frame
  } else {
    missing_taxa_phy <- input_process(missing_taxa) # process input if it is newick or phylo
    if (inherits(missing_taxa_phy, "phylo")) {
      phylo_check(phy = dated_tree, dated = TRUE)
      dtINmt <- dated_tree$tip.labels %in% missing_taxa$tip.labels
      mtINdt <- missing_taxa$tip.labels %in% dated_tree$tip.labels
      if (!all(dtINmt)) {
        warning("not all taxa from dated_tree are in missing_taxa tree")
      }
      missing_taxa_pruned <- ape::drop.tip(missing_taxa, missing_taxa$tip.labels[mtINdt])
      # phylo_prune_missing_taxa(phy = , taxa = ) # use this one??
      # dated_tree == missing_taxa_pruned # check that both trees are equal?
      # we don't need to
      # we can congruify if tree are not equal
      # we just need to make a tree with all lineages on dated_tree and all lineages in missing taxa vector
      missing_taxa <- missing_taxa_phy
      badformat <- FALSE
    } else {
      # if missing_taxa is provided as a vector, it must be a character vector:
      # it does not matter if original vector has no real names, if you want to add numbers or NAs or booleans at random, you can
      missing_taxa <- as.character(missing_taxa)
      missing_taxa[which(is.na(missing_taxa))] <- "NA"
      badformat <- FALSE
    }
  }
  if (length(missing_taxa) == 0) { # this is only valid if we don't wanna accept NULL as missing_taxa, but tere are cases in which it is, e.g. tree_fix_brlen
    missing_taxa <- NULL
    badformat <- FALSE
  }
  if (badformat) {
    stop("missing_taxa must be a character vector with species names,
		a data frame with taxonomic assignations, a newick character string, a phylo object, or NULL")
  }
  # IMPORTANT: Add a check that taxa in dated.trees is in reference_tree and viceversa
  return(missing_taxa)
  # if (is.null(reference_tree)){
  # if (is.null(taxon_summary)) {
  # cat("specify a reference_tree or taxon_summary to be added to the dated.trees")
  # stop("")
  # } else {
  # # construct a tree with phylotastic or take that from OpenTree?
  # cat("Constructing a reference_tree with taxa from dated.tree and taxon_summary")
  # }
  # } else {
  # if(!is.null(taxon_summary)){
  # cat("A reference_tree was given, taxon_summary argument is ignored")
  # }
}
