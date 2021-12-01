#' Date a given tree topology using a given set of calibrations
#'
#' @description \code{use_all_calibrations} generates one or multiple chronograms (i.e., phylogenetic
#' trees with branch lengths proportional to time) by dating a tree topology given
#' in \code{phy}, and secondary calibrations given in \code{calibrations},
#' using bladj or PATHd8.
#' @inheritParams use_calibrations_bladj
#' @param calibrations An object of class \code{datelifeCalibrations},
#' typically an output of [get_all_calibrations()].
#' @inheritParams extract_calibrations_phylo
#' @param dating_method Tree dating method. Options are "bladj" or "pathd8"
#' @inherit datelife_use return
#' @inherit use_calibrations details
#' @export
use_all_calibrations <- function(phy = NULL,
                                 calibrations = NULL,
                                 each = FALSE,
                                 dating_method = "bladj") {
    #
    # date the tree with bladj, pathd8 if branch lengths:
    # TODO
    # find a way to inherit params for use_calibrations_pathd8 and use_calibrations_bladj
    # phytools method for argument assignation should work, for example:
    # if (methods::hasArg("expand"))
    # if (methods::hasArg("giveup"))
    # if (methods::hasArg("type"))
    # if (methods::hasArg("root_age"))
    # if (methods::hasArg("match_calibrations"))
    if (each) {  # if calibrations is a list of data frames:
      res <- use_calibrations_each(phy = phy,
                                   calibrations = calibrations,
                                   dating_method = dating_method)
    } else {
      res <- use_calibrations(phy,
                              calibrations,
                              dating_method = dating_method)
    }
    # TODO: make a 'datelife' class for chronograms with calibrations data.frame attached, produced by datelife
    return(res)
}
#' Date a given tree topology by using a given list of calibrations independently,
#' to generate multiple hypothesis of time of divergence
#'
#' @description \code{use_calibrations_each} wraps \code{use_calibrations} to take each set of
#' given calibrations and use it independently as constraints for BLADJ or PATHd8
#' to date a given tree topology.
#' @inheritParams use_all_calibrations
#' @param ... Arguments to pass to \code{use_calibrations}.
#' @return A \code{multiPhylo} object of trees with branch lengths proportional
#' to time.
#' @inherit use_calibrations details
#' @export
# You can get the datelifeResult object or the list of chronograms first
# phylo_all1 <- get_datelife_result(input = my_phy)
# phylo_all2 <- summarize_datelife_result(phylo_all1)
# Either will work the same:
# use_calibrations_each(phy = my_phy, phylo_all = phylo_all1)
# use_calibrations_each(phy = my_phy, phylo_all = phylo_all2)
# You can also get the list of calibrations outside
# my_calibrations <- get_all_calibrations(input = phylo_all1, each = TRUE)
# use_calibrations_each(phy = my_phy, calibrations = my_calibrations)
# All this means that you can use your own set of calibrations, not necessarily the ones found only in datelife.
# phy <- tax_otolall[[1]]
# calibrations <- tax_eachcalall[[1]]
use_calibrations_each <- function(phy = NULL,
                                 calibrations = NULL,
                                 ...) {

		# date the tree with bladj, or pathd8 if branch lengths:
    if (inherits(calibrations, "data.frame")) {
      calibrations <- structure(list(calibrations), class = c("list", "datelifeCalibrations"))
    }
		chronograms <- lapply(calibrations, function(x)
                                        use_calibrations(phy = phy,
                                                         calibrations = x,
                                                         ...))
		class(chronograms) <- c("multiPhylo", "datelife")
		names(chronograms) <- names(calibrations)
		return(chronograms = chronograms)
}
#' Date a given tree topology using a combined set of given calibrations
#'
#' @description \code{use_calibrations} combines all given calibrations and uses them as
#' constraints to perform a dating analysis on a given tree topology, using BLADJ
#' if it has no branch lengths, or PATHd8 if the given tree topology has initial
#' branch lengths.
#'
#' @inheritParams use_all_calibrations
#' @inheritDotParams use_calibrations_pathd8
#' @return A \code{phylo} object with branch lengths proportional to time.
#' @export
#' @details
#' \code{calibrations} and \code{dating method} are stored as attributes.
#' They can be accessed with \code{attributes(return)$datelife_calibrations} and
#' \code{attributes(return)$dating_method}, respectively. \cr
#' If \code{input} tree does not have branch lengths, \code{dating_method} is
#' ignored and the function will use secondary calibrations to date the
#' \code{input} tree with [BLADJ](http://phylodiversity.net/bladj/),
#' using [phylocomr::ph_bladj()], via [use_calibrations_bladj()]. If
#' \code{input} tree has have branch lengths it can be dated
#' with [PATHd8](https://www2.math.su.se/PATHd8/), using [geiger::PATHd8.phylo()],
#' via [use_calibrations_pathd8()].
use_calibrations <- function(phy = NULL,
                             calibrations = NULL,
                             dating_method = "bladj",
                             type = "median",
                             ...){
	# check that phy names are in calibrations.df: done in match_all_calibrations inside use_calibrations_bladj
  message("... Using calibrations to date a tree topology.")
  exit <- FALSE
  if (inherits(phy, "multiPhylo")) {
    message(message_multiphylo())
    phy <- phy[[1]]
  }
	if(!inherits(phy, "phylo")){
    exit <- TRUE
		msg1 <- "'phy' is NOT a phylo object. Check this."
	} else {
    msg1 <- "'phy' is a phylo object. You are good."
  }
	# enhance: add a check for calibrations object structure
	if(!inherits(calibrations, "data.frame")){
    exit <- TRUE
		msg2 <- "'calibrations' is NOT a data frame. Check this."
	} else {
    msg2 <- "'calibrations' is a data frame. You are good."
  }
  if(exit){
    stop("This function requires valid values for both arguments 'phy' and 'calibrations'.\n",
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
  attr(chronogram, "datelife_calibrations") <- calibrations
  attr(chronogram, "dating_method") <- dating_method
  class(chronogram) <- c(class(chronogram), "datelife")
  message("Success!")
	return(chronogram)
}
