#' Perform a dating analysis on a set of input taxa.
#' It uses all calibrations from chronograms in the datelife database.
#' @inheritParams datelife_search
#' @inheritDotParams get_all_calibrations
#' @return A list with a chronogram and all calibrations found
#' @export
#' @details
#' If input is a tree, it will get all calibrations and date it with bladj if it
#' has no branch lengths or with PATHd8 if it has branch lengths.
#' If input is a vector of names, it will try to construct a bold tree first to get
#' a topology with branch lengths. If it can't, it will get the open tree of life
#' topology only and date it with bladj.
use_all_calibrations <- function(input = NULL, ...) { # dating_method = "bladj",
		# use congruification to expand calibrations: already implemented in map_all_calibrations
		# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]
		# this is just to be able to run if input is NULL, as an example:
		if(is.null(input)){
			phy <- make_bold_otol_tree(c("Rhea americana", "Struthio camelus", "Gallus gallus"), chronogram = FALSE, verbose = FALSE)
			if(!inherits(phy, "phylo")){
				phy <- get_dated_otol_induced_subtree(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
			}
		}
		datelife_query <- make_datelife_query(input = input, verbose = verbose) # it removes singleton nodes in phy
		# if input is not a tree, get one with bold or just otol:
		if(!inherits(datelife_query$phy, "phylo")){
			phy <- make_bold_otol_tree(datelife_query$cleaned_names, chronogram = FALSE, verbose = FALSE)
		} else {
			phy <- datelife_query$phy
		}
		# perform the datelife_search through get_all_calibrations:
		calibrations <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), ...)
		# date the tree with bladj, pathd8 if branch lengths:
		chronogram <- use_calibrations(phy, calibrations)
		return(list(phy = chronogram, calibrations.df = calibrations.df))
}

#' Perform a dating analysis on a tree topology using a determined set of calibrations.
#' @inheritParams use_calibrations_bladj
#' @param dating_method The method used for tree dating.
#' @inheritDotParams use_calibrations_pathd8
#' @return A phylo object with branch lengths propotional to time.
#' @export
#' @details
#' If phy does not have branch lengths, dating_method is ignored and BLADJ will be used.
use_calibrations <- function(phy = NULL, calibrations = NULL, dating_method = "bladj", type = "median", ...){
	# check that input names are in calibrations.df: done in map_all_calibrations inside use_calibrations_bladj
	if(!inherits(phy, "phylo")){
		message("phy is not a phylo object")
		return(NA)
	}
	# enhance: add a check for calibrations object structure
	dating_method <- match.arg(tolower(dating_method), c("bladj", "pathd8"))
	if(dating_method %in% "bladj"){
		phy$edge.length <- NULL
		chronogram <- use_all_calibrations_bladj(phy, calibrations, type)
	}
	if(dating_method %in% "pathd8"){
		if(is.null(phy$edge.length)){
			message('phy has no branch lengths, using dating_method = "bladj" instead.')
			chronogram <- use_all_calibrations_bladj(phy, calibrations, type)
		}
		chronogram <- use_all_calibrations_pathd8(phy, calibrations, ...)
	}
	return(chronogram)
}
