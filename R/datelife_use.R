#' Generate one or multiple chronograms for a set of given taxa.
#'
#' \code{datelife_use} generates one or multiple chronograms (i.e., phylogenetic
#' trees with branch lengths proportional to time) for a set of \code{input} taxa,
#' dated with bladj or PATHd8, using secondary calibrations available for any pair
#' of \code{input} taxa, mined from the code{\link[=opentree_chronograms]{opentree_chronograms}}
#' local DateLife database.
#'
#' @aliases use_all_calibrations
#' @inheritParams datelife_search input
#' @inheritParams use_all_calibrations dating_method
#' @inheritDotParams get_all_calibrations
#' @return A list with a chronogram and secondary calibrations used for dating.
#' @export
#' @details
#' If input is a tree, it will use secondary calibrations to (1) date the tree with bladj
#' if it has no branch lengths, or (2) date the tree with PATHd8 if it has branch lengths.
#' If input is a vector of taxon names, it will attempt to reconstruct a BOLD tree first
#' to get a topology with branch lengths. If it can't, it will get an Open Tree
#' of Life synthetic tree topology and will date it with bladj.
datelife_use <- function(input = NULL,
                         dating_method = "bladj",
                         ...) {
		# use congruification to expand calibrations: already implemented in match_all_calibrations
		# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]

    # get_all_calibrations arguments
    if (methods::hasArg("each")) {
      each <- list(...)$each
      if(!is.logical(each)){
        stop("'each' argument must be one of TRUE or FALSE")
      }
    } else {
      each <- FALSE
    }
		# run with an example of three birds when input is NULL
		if(is.null(input)){
      stop("'input' is NULL")
		}
		# run with taxa provided by user
	  # if(inherits(input, "datelifeResult") | any(is.numeric(input))){
	  #   # a datelifeResult object is a list of chronograms from OpenTree matching at least 2 queried taxa
	  #   stop("'input' must be any of the following: \n\t 1) a character vector of taxon names, \n\t 2) a tree as 'phylo' object or newick character string, or \n\t 3) a 'datelifeQuery' object, see 'make_datelife_query' function.")
	  # }
  	if(inherits(input, "multiPhylo")) {
  		input <- input[[1]]
  		message("'input' is a multiPhylo object. Only the first 'phylo' object will be used.")
  	}
    datelife_query <- input

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

    # perform a datelife_search through get_all_calibrations:
    calibrations <- get_all_calibrations(input = input, each = each)

		res <- use_all_calibrations(phy = phy,
                                calibrations = calibrations,
                                dating_method = dating_method)
    return(res)
}

#' Generate one or multiple chronograms for a set of given taxon names.
#'
#' \code{datelife_use_vector} generates one or multiple chronograms (i.e., phylogenetic
#' trees with branch lengths proportional to time) for a set of \code{input} taxa,
#' dated with bladj or PATHd8, using secondary calibrations available for any pair
#' of \code{input} taxa, mined from the code{\link[=opentree_chronograms]{opentree_chronograms}}
#' local DateLife database.
#'
#' @aliases use_all_calibrations_vector
#' @param input A character vector of taxon names
#' @inheritParams datelife_search
#' @param dating_method The method used for tree dating. Options are "bladj" or "pathd8"
#' @inheritDotParams get_all_calibrations
#' @return A list with a chronogram and secondary calibrations used
#' @export
#' @details
#' If input is a tree, it will use secondary calibrations to (1) date the tree with bladj
#' if it has no branch lengths, or (2) date the tree with PATHd8 if it has branch lengths.
#' If input is a vector of taxon names, it will attempt to reconstruct a BOLD tree first
#' to get a topology with branch lengths. If it can't, it will get an Open Tree
#' of Life synthetic tree topology and will date it with bladj.
datelife_use_vector <- function(input = NULL,
                         dating_method = "bladj",
                         ...) {

#
#
}
