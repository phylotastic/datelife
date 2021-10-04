#' Get a chronogram (dated phylogenetic tree) or published times of lineage divergence available for a set of taxa.

#' Core datelife function to input a vector of species, newick string, or phylo object to

#' @aliases datelife
#' @param input A character vector of taxon names, a tree as a 'phylo' object or
#' a newick character string, or a 'datelifeQuery' object from make_datelife_query function.
#' @param summary_format The desired output format for target chronograms (chronograms
#' of target taxa). See details.
#' @param summary_print A character vector specifying type of summary information
#' to be printed: "citations" for the references of chronograms from cache where target
#' taxa are found, "taxa" for a summary of the number of chronograms where each target
#' taxon is found, or "none" if nothing should be printed. Default to display both
#' c("citations", "taxa").
#' @param taxon_summary A character vector specifying if data on target taxa missing
#' in source chronograms should be added to the output as a "summary" or as a
#' presence/absence "matrix". Default to "none", no information on taxon_summary
#' added to the output.
#' @param partial Boolean; default to TRUE: use source trees even if they only
#' match some of the desired taxa.
#' @param use_tnrs Boolean; default to FALSE. If TRUE, use OpenTree's services
#' to resolve names. This can dramatically improve the chance of matches, but also
#' take longer.
#' @param approximate_match Boolean; default to TRUE: use a slower TNRS to correct
#' misspellings, increasing the chance of matches (including false matches)
#' @param update_opentree_chronograms default to FALSE
#' @param cache A character vector of length one, with the name of the data object
#' to cache. Default to "opentree_chronograms", a data object storing Open Tree of
#' Life's database chronograms and other associated information.
#' @param get_spp_from_taxon boolean vector, default to FALSE. If TRUE, will get
#' all species names from taxon names given in input. Must have same length as input.
#' If input is a newick string , with some clades it will be converted to phylo
#' object phy, and the order of get_spp_from_taxon will match phy$tip.label.
#' @param verbose Boolean. If TRUE, it gives printed updates to the user.
#' @param criterion Whether to get the grove with the most trees or the most taxa
#' @export
#' @details
#' Available output formats are:
#'
#' citations: A character vector of references where chronograms with some or all of the target taxa are published (source chronograms).
#'
#' mrca: A named numeric vector of most recent common ancestor (mrca) ages of target taxa defined in input, obtained from the source chronograms. Names of mrca vector are equal to citations.
#'
#' newick_all: A named character vector of newick strings corresponding to target chronograms derived from source chronograms. Names of newick_all vector are equal to citations.
#'
#' newick_sdm: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' newick_median: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram from the median of all source chronograms.
#'
#' phylo_sdm: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' phylo_median: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained from source chronograms with median method.
#'
#' phylo_all: A named list of phylo objects corresponding to each target chronogram obtained from available source chronograms. Names of phylo_all list correspond to citations.
#'
#' phylo_biggest: The chronogram with the most taxa. In the case of a tie, the chronogram with clade age closest to the median age of the equally large trees is returned.
#'
#' html: A character vector with an html string that can be saved and then opened in any web browser. It contains a 4 column table with data on target taxa: mrca, number of taxa, citations of source chronogram and newick target chronogram.
#'
#' data_frame A data frame with data on target taxa: mrca, number of taxa, citations of source chronograms and newick string.
#'
#' For approaches that return a single synthetic tree, it is important that the trees leading to it form a grove (roughly, a sufficiently overlapping set of taxa between trees: see Ané et al. 2005, 10.1007/s00026-009-0017-x). In the rare case of multiple groves, should we take the one with the most trees or the most taxa?
#'
#' @examples
#' # obtain median ages from a set of source chronograms in newick format:
#' ages <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="newick_median")
#' # save the tree in newick format
#' write(ages, file="some.bird.ages.txt")
#'
#' # obtain median ages from a set of source chronograms in phylo format
#' # will produce same tree as above but in r phylo format:
#' ages.again <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="phylo_median")
#' plot(ages.again)
#' library(ape)
#' ape::axisPhylo()
#' mtext("Time (million years ago)", side = 1, line = 2, at = (max(get("last_plot.phylo",
#' 		envir = .PlotPhyloEnv)$xx) * 0.5))
#' write.tree(ages.again, file="some.bird.tree.again.txt") # saves phylo object in newick format
#'
#' # obtain mrca ages and target chronograms from all source chronograms
#' # generate an html  output readable in any web browser:
#' ages.html <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="html")
#' write(ages.html, file="some.bird.trees.html")
#' system("open some.bird.trees.html")
datelife_search <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
														summary_format = "phylo_all",
														partial = TRUE,
														use_tnrs = FALSE,
														approximate_match = TRUE,
														update_opentree_chronograms = FALSE,
														cache = "opentree_chronograms",
														summary_print= c("citations", "taxa"),
														taxon_summary = c("none", "summary", "matrix"),
														get_spp_from_taxon = FALSE,
														verbose = FALSE,
														criterion="taxa") {
	# TODO: find a way not to repeat partial and cache arguments, which are used in both get_datelife_result and summarize_datelife_result
	if(update_opentree_chronograms){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	} else {
		if("opentree_chronograms" %in% cache){
			utils::data("opentree_chronograms")
			cache <- get("opentree_chronograms")
		}
	}
	datelife_query <- make_datelife_query(input = input, use_tnrs = use_tnrs,
		approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon,
		verbose = verbose)
	datelife_result.here <- get_datelife_result(input = datelife_query,
		partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match,
		update_opentree_chronograms = FALSE, cache = cache, verbose = verbose)
	# print.datelife(datelife_result = datelife_result.here)
	# datelife <- list(datelife_query = datelife_query, datelife_result = datelife_result.here)
	# class(datelife) <- "datelife"
	# return(datelife)
	return(summarize_datelife_result(datelife_result = datelife_result.here,
		datelife_query = datelife_query, summary_format = summary_format,
		partial = partial, update_opentree_chronograms = FALSE, cache = cache,
		summary_print = summary_print, taxon_summary = taxon_summary,
		verbose = verbose, criterion = criterion))
}

# print.datelife <- function(datelife){
# 	datelife_result_check(datelife$datelife_result)
# 	cat("Number and which queried taxa found in chronogram database")
# 	cat("Number of queried taxa not found")
# 	cat("Number of chronograms with at least two queried taxa found in database")
# }

#' Go from a vector of species, newick string, or phylo object to a list of patristic matrices
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
# #' @inheritParams make_bold_otol_tree
# #' @inheritDotParams make_bold_otol_tree
#' @return A datelifeResult object - a named list of patristic matrices. You can ccess the original query with attributes(my_datelife_result)$query
#' @export
get_datelife_result <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
																partial = TRUE,
																use_tnrs = FALSE,
																approximate_match = TRUE,
																update_opentree_chronograms = FALSE,
																cache = get("opentree_chronograms"),
																get_spp_from_taxon = FALSE,
																verbose = FALSE) {
	if(update_opentree_chronograms){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	if(is_datelife_query(input)){
		input_dq <- input
	} else {
		input_dq <- make_datelife_query(input = input, use_tnrs = use_tnrs,
			approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon,
			verbose = verbose)
	}
	if(length(input_dq$cleaned_names) == 1){
			message("Cannot perform a search of divergence times with just one taxon.")
			if(!get_spp_from_taxon) {
				# message("Performing a clade search?? set get_spp_from_taxon = TRUE")
				message("Setting up get_spp_from_taxon = TRUE")
				input_dq <- make_datelife_query(input = input_dq$cleaned_names,
					get_spp_from_taxon = TRUE, use_tnrs = use_tnrs,
					approximate_match = approximate_match, verbose = verbose)
			}
	}
	if(length(input_dq$cleaned_names) == 1){
		message("Input has length one (even after searching spp within clades).")
		warning("\t", "Input contains only one lineage.")
		return(NA)
	}
	# setting phy to NULL always; it is a bad idea to congruidy subset trees,
	# do that later in summarizing steps
	results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = input_dq$cleaned_names, phy = NULL)
  	datelife_result <- results_list_process(results_list, input_dq$cleaned_names, partial)
	datelife_result_check(datelife_result, use_tnrs)
	class(datelife_result) <- c("datelifeResult")
	attr(datelife_result, "query") <- input_dq
	return(datelife_result)
}

#' checks if we obtained an empty search with the set of input taxon names
#' @inheritParams datelife_search
#' @param datelife_result An output of get_datelife_result function: A list of patristic matrices with names corresponding to original study citations.
#' @export
datelife_result_check <- function(datelife_result, use_tnrs = FALSE, verbose = FALSE){
	if(length(datelife_result) < 1) {
		warning("Datelife Result object is empty.", call. = FALSE)
		# if(verbose) {
			message("Input taxa were not found across available chronograms.")
			if(!use_tnrs) {
				message("Setting use_tnrs = TRUE might change this, but it is time consuming.")
			}
		# }
	}
}
#' checks if we obtained an empty search with the set of input taxon names
#' @inheritParams datelife_result_check
#' @return a datelifeResult object or a message explaining why it is not a datelifeResult object
check_datelife_result <- function(datelife_result){
	# datelife_result <- subset2_result
	if(!inherits(datelife_result, "datelifeResult")){
		if(is.list(datelife_result)){
			class_dl <- unname(sapply(datelife_result, class))
			if(all(grepl("matrix", class_dl))){
				class(datelife_result) <- "datelifeResult"
			} else {
				message("datelife_result has some elements that are not matrices; check it out.")
			}
		} else {
			message("datelife_result is not a list; check it out")
		}
	}
	return(datelife_result)
}

#' datelifeResult object of three birds "Rhea americana", "Pterocnemia pennata", and "Struthio camelus"
#'
#' @name threebirds_dr
#' @docType data
#' @format A list of 9 named patristic matrix
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree subset chronogram ants datelife
#' @details
#' Generated with:
#' threebirds_dr <- get_datelife_result(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
#' partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, cache = "opentree_chronograms")
#' use_data(threebirds_dr)
"threebirds_dr"
