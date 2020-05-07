#' Perform a dating analysis on a set of input taxa.
#' It searches and uses all calibrations from chronograms in the DateLife database.
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
		# use congruification to expand calibrations: already implemented in match_all_calibrations
		# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]
		# this is just to be able to run if input is NULL, as an example:
		if(is.null(input)){ # run with an example of three birds
			phy <- make_bold_otol_tree(c("Rhea americana", "Struthio camelus", "Gallus gallus"), chronogram = FALSE, verbose = FALSE)
			if(!inherits(phy, "phylo")){
				phy <- get_dated_otol_induced_subtree(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
			}
		} else { # run with taxa provided by user
			datelife_query <- datelife_query_check(datelife_query = input) # this also removes singleton nodes in phy
			# if input is not a tree, get one with bold or just otol:
			if(!inherits(datelife_query$phy, "phylo")){
				phy <- make_bold_otol_tree(datelife_query$cleaned_names, chronogram = FALSE, verbose = FALSE)
			} else {
				phy <- datelife_query$phy
			}
		}

		# perform the datelife_search through get_all_calibrations:
		calibrations <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), ...)
		# date the tree with bladj, pathd8 if branch lengths:
		chronogram <- use_calibrations(phy, calibrations)
		return(list(phy = chronogram, calibrations.df = calibrations))
}

#' Use calibrations from each source chronogram to date a tree.
#' @param phy A tree with or without branch lengths to be dated. Only as phylo object for now.
#' @param phylo_all Can be NULL. A datelifeResult list of patristic matrices, or a chronogram as phylo or multiPhylo.
#' @param calibrations Can be NULL. A list of dataframes from get_all_calibrations function.
#' @param ... Arguments to pass to get_all_calibrations
#' @return A list with a multiPhylo object of chronograms and a list of all calibrations used
#' @export
#' @details
#' You can  get the datelifeResult object or the list of chronograms first
#' phylo_all1 <- get_datelife_result(input = my_phy)
#' phylo_all2 <- summarize_datelife_result(phylo_all1)
#' Either will work the same:
#' use_each_calibration(phy = my_phy, phylo_all = phylo_all1)
#' use_each_calibration(phy = my_phy, phylo_all = phylo_all2)
#' You can also get the list of calibrations outside
#' my_calibrations <- get_all_calibrations(input = phylo_all1, each = TRUE)
#' use_each_calibration(phy = my_phy, calibrations = my_calibrations)
#' All this means that you can use your own set of calibrations, not necessarily the ones found only in datelife.
# phy <- tax_otolall[[1]]
# calibrations <- tax_eachcalall[[1]]
use_each_calibration <- function(phy = NULL, phylo_all = NULL, calibrations = NULL, ...) {
		if(!inherits(phy, "phylo")){
			message("phy is not a phylo object")
			return(NA)
		}
		# when calibrations are not given as a list of dataframes
		if(!inherits(calibrations, "list")){
			if(inherits(phylo_all, "datelifeResult") | inherits(phylo_all, "phylo") | inherits(phylo_all, "multiPhylo")){
				calibrations <- get_all_calibrations(input = phylo_all, each = TRUE, ...)
			} else {
				# can perform the datelife_search through get_all_calibrations:
				calibrations <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), each = TRUE, ...)
			}
		}
		# date the tree with bladj, or pathd8 if branch lengths:
		chronograms <- lapply(calibrations, function(x) use_calibrations(phy, x))
		class(chronograms) <- "multiPhylo"
		names(chronograms) <- names(calibrations)
		return(list(phy = chronograms, calibrations = calibrations))
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
	# check that input names are in calibrations.df: done in match_all_calibrations inside use_calibrations_bladj
	if(!inherits(phy, "phylo")){
		message("phy is not a phylo object")
		return(NA)
	}
	# enhance: add a check for calibrations object structure
	dating_method <- match.arg(tolower(dating_method), c("bladj", "pathd8"))
	if(dating_method %in% "bladj"){
		phy$edge.length <- NULL
		chronogram <- use_calibrations_bladj(phy, calibrations, type)
	}
	if(dating_method %in% "pathd8"){
		if(is.null(phy$edge.length)){
			message('phy has no branch lengths, using dating_method = "bladj" instead.')
			chronogram <- use_calibrations_bladj(phy, calibrations, type)
		}
		chronogram <- use_calibrations_pathd8(phy, calibrations, ...)
	}
	return(chronogram)
}
