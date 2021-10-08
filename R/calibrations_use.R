#' Use given calibrations to generate one or multiple chronograms for a given
#' tree topology.
#'
#' \code{use_calibrations_each} uses each set of calibrations given by the user
#' independently as constraints for BLADJ or PATHd8 to date a given tree topology.
#'
#' @param phy A tree with or without branch lengths to be dated. Only as phylo object for now.
#' @param calibrations Can be NULL. A list of dataframes from get_all_calibrations function.
#' @param ... Arguments to pass to use_calibrations
#' @return A list with a multiPhylo object of chronograms and a list of all calibrations used
#' @export
#' @details
#' You can get the datelifeResult object or the list of chronograms first
#' phylo_all1 <- get_datelife_result(input = my_phy)
#' phylo_all2 <- summarize_datelife_result(phylo_all1)
#' Either will work the same:
#' use_calibrations_each(phy = my_phy, phylo_all = phylo_all1)
#' use_calibrations_each(phy = my_phy, phylo_all = phylo_all2)
#' You can also get the list of calibrations outside
#' my_calibrations <- get_all_calibrations(input = phylo_all1, each = TRUE)
#' use_calibrations_each(phy = my_phy, calibrations = my_calibrations)
#' All this means that you can use your own set of calibrations, not necessarily the ones found only in datelife.
# phy <- tax_otolall[[1]]
# calibrations <- tax_eachcalall[[1]]
use_calibrations_each <- function(phy = NULL,
                                 calibrations = NULL,
                                 ...) {

		# date the tree with bladj, or pathd8 if branch lengths:
		chronograms <- lapply(calibrations, function(x) use_calibrations(phy, x, ...))
		class(chronograms) <- "multiPhylo"
		names(chronograms) <- names(calibrations)
		return(list(chronograms = chronograms, calibrations = calibrations))
}
#' Use given calibrations to generate a chronogram for a given tree topology.
#'
#' \code{use_calibrations} combines all given calibrations given and uses them as
#' constraints to perform a dating analysis on a given tree topology, using BLADJ or PATHd8.
#'
#' @inheritParams use_calibrations_bladj
#' @inheritParams datelife
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
