# Functions to perform a molecular dating analysis on a set of input taxa, using calibrations from chronograms in a database.

#' Use all calibrations from chronograms in a database to date a tree.
#' @param phy A phylo object
#' @param all_calibrations A data frame of calibrations from get_all_calibrations function
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams datelife_search
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
use_all_calibrations <- function(phy = NULL,
	all_calibrations = NULL, partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE,
	update_cache = FALSE, cache = get("opentree_chronograms"), verbose = FALSE) {
		# enhance: use congruification to expand calibrations:
			# already implemented in map_all_calibrations(phy, get_all_caibrations(phy))$calibrations
			# and pathd8 still does not work sometimes
		# calibrations.df <- eachcal[[2]]
		# calibrations.df <- calibs$calibration
		# phy <- tax_phyloallall[[2]][[3]]
		if(is.null(phy)){ # just to run an example:
			phy <- make_bold_otol_tree(c("Rhea americana", "Struthio camelus", "Gallus gallus"), chronogram = TRUE, verbose = FALSE)
			if(!inherits(phy, "phylo")){
				phy <- get_dated_otol_induced_subtree(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
			}
		} else {
			input <- input_process(input = phy, verbose = verbose)
			phy <- input
		}
		# phy must be a tree, check this
		if(!inherits(phy, "phylo")){
			message("phy is not a phylo object")
			return(NA)
		}
		# remove singleton nodes in phy:
		phy <- ape::collapse.singles(phy)
		if(is.null(all_calibrations)){
			calibrations.df <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label),
			partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match,
			update_cache = update_cache, cache = cache, verbose = verbose)
		} else {
			# enhance: add a check for object structure
			# enhance: check that input names are in calibrations.df
			calibrations.df <- all_calibrations
		}

		return(list(phy = chronogram, calibrations.df = calibrations.df))
}
