#' Generate one or multiple chronograms for a set of given taxa
#'
#' @description \code{datelife_use} generates one or multiple chronograms (i.e., phylogenetic
#' trees with branch lengths proportional to time) for a set of \code{input} taxa,
#' dated with bladj or PATHd8, using secondary calibrations available for any pair
#' of \code{input} taxa, mined from the code{\link[=opentree_chronograms]{opentree_chronograms}}
#' local DateLife database.
#'
#' @aliases use_all_calibrations
#' @inheritParams datelife_search
#' @inheritParams use_all_calibrations
#' @return A list with a chronogram and secondary calibrations used for dating.
#' @export
#' @details
#' If input is a tree, it will use secondary calibrations to (1) date the tree with bladj
#' if it has no branch lengths, or (2) date the tree with PATHd8 if it has branch lengths.
#' If input is a vector of taxon names, it will attempt to reconstruct a BOLD tree first
#' to get a topology with branch lengths. If it can't, it will get an Open Tree
#' of Life synthetic tree topology and will date it with bladj.
datelife_use <- function(input = NULL,
                         each = FALSE,
                         dating_method = "bladj") {
		# use congruification to expand calibrations: already implemented in match_all_calibrations
		# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]

		# run with taxa provided by user
	  # if(inherits(input, "datelifeResult") | any(is.numeric(input))){
	  #   # a datelifeResult object is a list of chronograms from OpenTree matching at least 2 queried taxa
	  #   stop("'input' must be any of the following: \n\t 1) a character vector of taxon names, \n\t 2) a tree as 'phylo' object or newick character string, or \n\t 3) a 'datelifeQuery' object, see 'make_datelife_query' function.")
	  # }


    if(!is_datelife_query(datelife_query)){
      # make_datelife_query also removes singleton nodes in phy
      # TODO: add extra arguments for make_datelife_query function with hasArg (phytools method)
      datelife_query <- make_datelife_query(input = input)
    }
	  # if input is not a tree, get one with bold or otol:
	  if(!inherits(datelife_query$phy, "phylo")){
	  	phy <- make_bold_otol_tree(datelife_query$cleaned_names, chronogram = FALSE, verbose = FALSE)
	  	if(!inherits(phy, "phylo")){
        message("BOLD tree reconstruction failed for this set of taxa.")
	  		phy <- get_dated_otol_induced_subtree(input = datelife_query$cleaned_names)
	  	}
	  } else {
	  	phy <- datelife_query$phy
	  }

    # call datelife_search, get_all_calibrations and use_all_calibrations:
    res <- datelife_use_phylo(phy = phy,
                              each = each,
                              dating_method = dating_method)

    return(res)
}

#' Generate one or multiple chronograms for a set of given taxon names
#'
#' \code{datelife_use_basic} generates one or multiple chronograms (i.e., phylogenetic
#' trees with branch lengths proportional to time) for a set of \code{input} taxa,
#' dated with bladj or PATHd8, using secondary calibrations available for any pair
#' of \code{input} taxa, mined from the code{\link[=opentree_chronograms]{opentree_chronograms}}
#' local DateLife database.
#'
#' @param input A \code{datelifeQuery} object.
#' @inheritDotParams get_all_calibrations
#' @inheritDotParams use_all_calibrations
#' @return A list with a chronogram and secondary calibrations used.
#' @export
datelife_use_datelifequery <- function(input = NULL,
                                       dating_method = "bladj",
                                       each = FALSE) {
#
    if (!is_datelife_query(input)) {
      warning("'input' is not a 'datelifeQuery' object.")
      return(NA)
    }
    if (!inherits(input$phy, phylo)) {
      warning("'input$phy' is not a 'phylo' object.")
      return(NA)
    }
    # get calibrations by performing a datelife_search with get_all_calibrations:
    calibrations <- get_calibrations_datelifequery(input = input,
                                                   each = each)

    # date the topology with obtained calibrations
    res <- use_all_calibrations(phy = input$phy,
                                calibrations = calibrations,
                                dating_method = dating_method)
    return(res)
}
#'
#' @param phy a \code{phylo} object with or without branch lengths.
#' @inheritParams get_all_calibrations
#' @inherit datelife_use return
#' @export
datelife_use_phylo <- function(phy = NULL,
                               each = FALSE,
                               dating_method = "bladj") {
    #
    if(inherits(phy, "multiPhylo")) {
      message(message_multiphylo())
      phy <- phy[[1]]
    }
    if(!inherits(phy, "phylo")) {
      warning("'phy' must be a 'phylo' object.")
      return(NA)
    }
    # get_all_calibrations uses phy$tip.label to do a datelife_search, make a
    # datelifeQuery, and extract calibrations:
    calibrations <- get_all_calibrations(input = phy,
                                         each = each)

    # date the given tree topology with calibrations obtained in previous step
    res <- use_all_calibrations(phy = phy,
                                calibrations = calibrations,
                                dating_method = dating_method)
    return(res)
}
