# Functions to perform a molecular dating analysis on a set of input taxa, using calibrations from chronograms in a database.

#' Use all calibrations from chronograms in a database to date a tree.
#' @param phy A phylo object
#' @param all_calibrations A data frame of calibrations from get_all_calibrations function
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param expand How much to expand by each step to get consistent calibrations
#' @param giveup How many expansions to try before giving up
#' @inheritParams datelife_search
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
use_all_calibrations <- function(phy = NULL,
	all_calibrations = NULL, partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE,
	update_cache = FALSE, cache = get("opentree_chronograms"), expand = 0.1,
	giveup = 100, verbose = FALSE) {
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
		phy$tip.label <- gsub(' ', '_', phy$tip.label) #underscores vs spaces: the battle will never end.
		calibrations.df$taxonA <- gsub(' ', '_', calibrations.df$taxonA)
		calibrations.df$taxonB <- gsub(' ', '_', calibrations.df$taxonB)
		calibrations.df <- calibrations.df[which(calibrations.df$taxonA %in% phy$tip.label),]
		calibrations.df <- calibrations.df[which(calibrations.df$taxonB %in% phy$tip.label),]
		if(nrow(calibrations.df) == 0){
			message("phy cannot be dated")
			return(list(phy = NA, calibrations.df = calibrations.df))
		}
		chronogram <- NULL
		if(expand != 0) {
			attempts = 0
			original.calibrations.df <- calibrations.df
			try(chronogram <- geiger::PATHd8.phylo(phy, calibrations.df), silent = TRUE)
			# original.calibrations.df <- get_all_calibrations(phy2)
			while(is.null(chronogram) & attempts < giveup) {
				# print(rep(attempts, 10))
				calibrations.df <- original.calibrations.df
				calibrations.df$MaxAge <- calibrations.df$MaxAge * ((1+expand)^attempts)
				calibrations.df$MinAge <- calibrations.df$MinAge * ((1-expand)^attempts)
				# We will have no fixed ages. Pathd8 just quietly gives up. So instead, we add a tiny branch with a zero calibration
				# between it and its sister.
				made.up.edgelength <- min(1e-9, .001*min(phy$edge.length))
				phy2 <- phytools::bind.tip(ape::reorder.phylo(phy), "tinytip", edge.length = made.up.edgelength, where = 1, position = made.up.edgelength) #bind tip has weird behavior for non-reordered trees
				calibrations.df[dim(calibrations.df)[1]+1,]<- c("fixed", 0, 0, phy$tip.label[1], "tinytip", "none")
				chronogram <- tryCatch(geiger::PATHd8.phylo(phy2, calibrations.df), error = function(e) NULL)
				# data that errors in geiger::PATHd8.phylo(phy2, calibrations.df) (and hence the tryCatch):
				# phy2 <- cetaceae_phyloall[2]# cetaceae_phyloall[3] and cetaceae_phyloall[4] also error
				# calibrations.df <- get_all_calibrations(cetaceae_phyloall[2])
				# Error in write.pathd8(phy, calibrations, base): Some calibrations not encountered in tree
				if(!is.null(chronogram)) {
					chronogram$edge.length[which(chronogram$edge.length < 0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
					chronogram <- ape::drop.tip(chronogram, "tinytip")
				}
				attempts <- attempts+1
			}
			if(attempts > 0 & inherits(chronogram, "phylo")) {
				message("Dates are even more approximate than usual: had to expand constraints to have them agree")
			}
		} else { # alternative to expansion: summarize calibrations
			calibs <- map_all_calibrations(phy, calibrations.df)
			try(chronogram <- geiger::PATHd8.phylo(calibs$phy, calibs$calibrations), silent = TRUE)
			calibrations.df <- calibs$calibrations
		}
		if(!is.null(chronogram)) {
			# sometimes pathd8 returns tiny negative branch lengths.
			# https://github.com/phylotastic/datelife/issues/11
			chronogram <- tree_fix_brlen(tree = chronogram, fixing_criterion =
				"negative", fixing_method = 0, ultrametric = TRUE)
			# chronogram$edge.length[which(chronogram$edge.length<0)] <- 0
		}
		return(list(phy = chronogram, calibrations.df = calibrations.df))
}
