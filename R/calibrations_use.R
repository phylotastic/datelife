#' Get secondary calibrations to generate one or multiple chronograms for a set of input taxa.
#'
#' \code{use_all_calibrations} returns a phylogenetic tree with branch lengths
#' proportional to time (i.e., a chronogram) for all \code{input} taxa, dated with
#' bladj or PATHd8, using secondary calibrations available from a chronogram database.
#'
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
use_all_calibrations <- function(input = NULL,
                                 dating_method = "bladj",
                                 ...) {
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
		}
		# run with taxa provided by user
	  if(inherits(input, "datelifeResult") | any(is.numeric(input))){
	    # a datelifeResult object is a list of chronograms from OpenTree matching at least 2 queried taxa
	    stop("'input' must be any of the following: \n\t 1) a character vector of taxon names, \n\t 2) a tree as 'phylo' object or newick character string, or \n\t 3) a 'datelifeQuery' object, see 'make_datelife_query' function.")
	  }
  	if(inherits(input, "multiPhylo")) {
  		input <- input[[1]]
  		message("'input' is a multiPhylo object. Only the first 'phylo' object will be used.")
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

#' Get one or multiple chronograms for a set of input taxa.
#'
#' \code{use_each_calibration} uses calibrations from each source chronogram
#' independently to date a tree from a set of input taxa.
#'
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
use_each_calibration <- function(phy = NULL,
                                 phylo_all = NULL,
                                 calibrations = NULL,
                                 ...) {
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
#' Get a chronogram for a set of input taxa.
#'
#' \code{use_calibrations} combines all calibrations from source chronograms
#' to date a tree from a set of input taxa.
#'
#' Perform a dating analysis on a tree topology using a given set of calibrations.
#' @inheritParams use_calibrations_bladj
#' @inheritParams use_all_calibrations
#' @inheritDotParams use_calibrations_pathd8
#' @return A phylo object with branch lengths proportional to time.
#' @export
#' @details
#' If phy does not have branch lengths, dating_method is ignored and BLADJ will be used.
use_calibrations <- function(phy = NULL,
                             calibrations = NULL,
                             dating_method = "bladj",
                             type = "median",
                             ...){
	# check that input names are in calibrations.df: done in match_all_calibrations inside use_calibrations_bladj
  exit <- FALSE
	if(!inherits(phy, "phylo")){
    exit <- TRUE
		msg1 <- "Value provided in 'phy' is NOT a phylo object. Check this."
	} else {
    msg1 <- "Value provided in 'phy' is a phylo object. You are good."
  }
	# enhance: add a check for calibrations object structure
	if(!inherits(calibrations, "data.frame")){
    exit <- TRUE
		msg2 <- "Value provided in 'calibrations' is NOT a data frame. Check this."
	} else {
    msg2 <- "Value provided in 'calibrations' is a data frame. You are good."
  }
  if(exit){
    stop("This function requires valid values for both arguments `phy` and `calibrations`.\n",
         msg1, "\n",
         msg2, "\n")
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
