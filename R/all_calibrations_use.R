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
use_all_calibrations <- function(input = NULL, dating_method = "bladj", ...) { # dating_method = "bladj",
		# use congruification to expand calibrations: already implemented in match_all_calibrations
		# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]

    # get_all_calibrations arguments
    if(methods::hasArg("each")){
      each <- list(...)$each
      if(!is.logical(each)){
        stop("'each' argument must be one of TRUE or FALSE")
      }
    } else {
      each <- FALSE
    }

		# run with an example of three birds when input is NULL
		if(is.null(input)){
			phy <- make_bold_otol_tree(c("Rhea americana", "Struthio camelus", "Gallus gallus"), chronogram = FALSE, verbose = FALSE)
			if(!inherits(phy, "phylo")){
				phy <- get_dated_otol_induced_subtree(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
			}
		} else {
			# run with taxa provided by user
		  if(inherits(input, "datelifeResult") | any(is.numeric(input))){
		    # a datelifeResult object is a list of chronograms from OpenTree matching at least 2 queried taxa
		    stop("Input must be any of the following: \n\t 1) a character vector of taxon names, \n\t 2) a tree as 'phylo' object or newick character string, or \n\t 3) a 'datelifeQuery' object, see 'make_datelife_query' function.")
		  } else {
      	if(inherits(input, "multiPhylo")) {
      		input <- input[[1]]
      		message("input is a multiPhylo object. Only the first 'phylo' object will be used.")
      	}
		    if(!is_datelife_query(input)){
		      message("Running 'make_datelife_query'...")
          # make_datelife_query also removes singleton nodes in phy
		      # TODO: add extra arguments for make_datelife_query function with hasArg (phytools method)
		      datelife_query <- make_datelife_query(input = input)
		    } else {
		      datelife_query <- input
		    }
			  # if input is not a tree, get one with bold or otol:
			  if(!inherits(datelife_query$phy, "phylo")){
			  	phy <- make_bold_otol_tree(datelife_query$cleaned_names, chronogram = FALSE, verbose = FALSE)
			  	if(!inherits(phy, "phylo")){
			  		phy <- get_dated_otol_induced_subtree(input = datelife_query$cleaned_names)
			  	}
			  } else {
			  	phy <- datelife_query$phy
			  }
		  }
		}

		# perform a datelife_search through get_all_calibrations:
		calibrations <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), each = each)
		# date the tree with bladj, pathd8 if branch lengths:
		# TODO
		# find a way to inherit params for use_calibrations_pathd8
		# phytools method for argument assignation seems a nice option

		#if calibrations is a list of data frames:
		if(each){
  		chronogram <- use_each_calibration(phy=phy, calibrations=calibrations, dating_method = dating_method)$chronograms
		} else {
			chronogram <- use_calibrations(phy, calibrations, dating_method = dating_method)
		}
		# TODO: make a 'datelife' class for chronograms with calibrations data.frame attached, produced by use_all_calibrations
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
		# when calibrations are not given as a list of data frames
		if(!inherits(calibrations, "list")){
			if(inherits(phylo_all, "datelifeResult") | inherits(phylo_all, "phylo") | inherits(phylo_all, "multiPhylo")){
				calibrations <- get_all_calibrations(input = phylo_all, each = TRUE)
			} else {
				# can perform the datelife_search through get_all_calibrations:
				calibrations <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), each = TRUE)
			}
		}
		# date the tree with bladj, or pathd8 if branch lengths:
		chronograms <- lapply(calibrations, function(x) use_calibrations(phy, x, ...))
		class(chronograms) <- "multiPhylo"
		names(chronograms) <- names(calibrations)
		return(list(chronograms = chronograms, calibrations = calibrations))
}
#' Perform a dating analysis on a tree topology using a given set of calibrations.
#' @inheritParams use_calibrations_bladj
#' @param dating_method The method used for tree dating.
#' @inheritDotParams use_calibrations_pathd8
#' @return A phylo object with branch lengths proportional to time.
#' @export
#' @details
#' If phy does not have branch lengths, dating_method is ignored and BLADJ will be used.
use_calibrations <- function(phy = NULL, calibrations = NULL, dating_method = "bladj", type = "median", ...){
	# check that input names are in calibrations.df: done in match_all_calibrations inside use_calibrations_bladj
	if(!inherits(phy, "phylo")){
		stop("'phy' is not a phylo object.")
	}
	# enhance: add a check for calibrations object structure
	if(!inherits(calibrations, "data.frame")){
		stop("'calibrations' is not a data.frame.\n\t Provide a set of calibrations for 'phy', hint: see get_all_calibrations function.")
	}
	dating_method <- match.arg(tolower(dating_method), c("bladj", "pathd8"))
	if(dating_method %in% "bladj"){
		phy$edge.length <- NULL
		chronogram <- use_calibrations_bladj(phy, calibrations, type)
	}
	if(dating_method %in% "pathd8"){
		if(is.null(phy$edge.length)){
			message("\nPATHd8 requires initial branch lengths and 'phy' has no branch lengths, using dating_method = bladj instead.\n")
			chronogram <- use_calibrations_bladj(phy, calibrations, type)
			dating_method <- "bladj"
		} else {
			chronogram <- use_calibrations_pathd8(phy, calibrations, ...)
		}
	}
	# TODO
	# add dating_method attribute to chronogram
	return(chronogram)
}
