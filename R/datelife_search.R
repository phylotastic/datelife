#' Get scientific, peer-reviewed information on time of lineage
#' divergence openly available for a given set of taxon names
#'
#' @description \code{datelife_search} is the core DateLife function to find and
#' get all openly available, peer-reviewed scientific information on time of
#' lineage divergence for a set \code{input} taxon names given as a character
#' vector, a newick character string, or a \code{phylo} or \code{multiPhylo} object.
#'
#' @aliases datelife
#' @param input One of the following:
#' \describe{
#'	\item{Taxon names}{As a character vector.}
#'	\item{One or more trees with taxon names as tip labels}{As a \code{phylo} or \code{multiPhylo}
#'				object, OR as a newick character string.}
#'	\item{A \code{datelifeQuery} object}{An output from [make_datelife_query()]
#' 				\code{\link[=make_datelife_query]{make_datelife_query}}.}
#'	}
#' @param summary_format A character vector of length one, indicating the output
#' format for results of the DateLife search. Available output formats are:
#' \describe{
#'	 \item{"citations"}{A character vector of references where chronograms with
#'	 				some or all of the target taxa are published (source chronograms).}
#'	 \item{"mrca"}{A named numeric vector of most recent common ancestor (mrca)
#'	 				ages of target taxa defined in input, obtained from the source chronograms.
#'	 				Names of mrca vector are equal to citations.}
#'	 \item{"newick_all"}{A named character vector of newick strings corresponding
#'	 				to target chronograms derived from source chronograms. Names of newick_all
#'	 				vector are equal to citations.}
#'	 \item{"newick_sdm"}{Only if multiple source chronograms are available. A
#'	 				character vector with a single newick string corresponding to a target
#'	 				chronogram obtained with SDM supertree method (Criscuolo et al. 2006).}
#'	 \item{"newick_median"}{Only if multiple source chronograms are available.
#'	 				A character vector with a single newick string corresponding to a target
#'	 				chronogram from the median of all source chronograms.}
#'	 \item{"phylo_sdm"}{Only if multiple source chronograms are available. A
#'	 				phylo object with a single target chronogram obtained with SDM supertree
#'	 				method (Criscuolo et al. 2006).}
#'	 \item{"phylo_median"}{Only if multiple source chronograms are available. A
#'	 				phylo object with a single target chronogram obtained from source
#'	 				chronograms with median method.}
#'	 \item{"phylo_all"}{A named list of phylo objects corresponding to each target
#'	 				chronogram obtained from available source chronograms. Names of
#'	 				phylo_all list correspond to citations.}
#'	 \item{"phylo_biggest"}{The chronogram with the most taxa. In the case of a
#'	 				tie, the chronogram with clade age closest to the median age of the
#'	 				equally large trees is returned.}
#'	 \item{"html"}{A character vector with an html string that can be saved and
#'	 				then opened in any web browser. It contains a 4 column table with data on
#'	 				target taxa: mrca, number of taxa, citations of source chronogram and
#'	 				newick target chronogram.}
#'	 \item{"data_frame"}{A data frame with data on target taxa: mrca, number of
#'	 				taxa, citations of source chronograms and newick string.}
#' }
#' @param summary_print A character vector specifying type of summary information
#' to be printed to screen:
#' \describe{
#'	 \item{"citations"}{Prints references of chronograms where target taxa are found.}
#'	 \item{"taxa"}{Prints a summary of the number of chronograms where each target
#'	 				taxon is found.}
#'	 \item{"none"}{Nothing is printed to screen.}
#' }
#' Default to \code{c("citations", "taxa")}, which displays both.
#' @param taxon_summary A character vector specifying if data on target taxa missing
#' in source chronograms should be added to the output as a \code{"summary"} or as a
#' presence/absence \code{"matrix"}. Default to \code{"none"}, no information on taxon_summary
#' added to the output.
#' @param partial Boolean; default to \code{TRUE}: use source trees even if they only
#' match some of the desired taxa.
#' @param update_opentree_chronograms default to \code{FALSE}
#' @param cache A character vector of length one, with the name of the data object
#' to cache. Default to \code{"opentree_chronograms"}, a data object storing Open Tree of
#' Life's database chronograms and other associated information.
#' @param verbose Boolean. If \code{TRUE}, it gives printed updates to the user.
#' @param criterion Whether to get the grove with the most trees or the most taxa
#' @export
#' @details
#' If only one taxon names is given as \code{input}, \code{get_spp_from_taxon} is
#' always set to \code{TRUE}.
#' For approaches that return a single synthetic tree, it is important that the
#' trees leading to it form a grove (roughly, a sufficiently overlapping set of
#' taxa between trees: see Ané et al. 2005, 10.1007/s00026-009-0017-x). In the
#' rare case of multiple groves, should we take the one with the most trees or
#' the most taxa?
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
														cache = "opentree_chronograms",
														summary_print= c("citations", "taxa"),
														taxon_summary = c("none", "summary", "matrix"),
														get_spp_from_taxon = FALSE,
														verbose = FALSE,
														criterion="taxa") {
	# TODO: find a way not to repeat partial and cache arguments, they are used in
	# datelife_search,  get_datelife_result and summarize_datelife_result
	# maybe take out option to update cache and work with versions of cache for reproducibility
	if("opentree_chronograms" %in% cache){
		utils::data("opentree_chronograms", package = "datelife")
		cache <- get("opentree_chronograms")
	}
	message("... Running a DateLife search.")
	datelife_query <- input
	if (suppressMessages(!is_datelife_query(input))){
		datelife_query <- make_datelife_query(input = input,
																					use_tnrs = use_tnrs,
																					approximate_match = approximate_match,
																					get_spp_from_taxon = get_spp_from_taxon,
																					verbose = verbose)
	}
	datelife_result.here <- get_datelife_result_datelifequery(input = datelife_query,
																														partial = partial,
																														use_tnrs = use_tnrs,
																														approximate_match = approximate_match,
																														cache = cache,
																														verbose = verbose)
	# print.datelife(datelife_result = datelife_result.here)
	# datelife <- list(datelife_query = datelife_query, datelife_result = datelife_result.here)
	# class(datelife) <- "datelife"
	# return(datelife)
	res <- summarize_datelife_result(datelife_result = datelife_result.here,
																	 datelife_query = datelife_query,
																	 summary_format = summary_format,
																	 partial = partial,
																	 cache = cache,
																	 summary_print = summary_print,
																	 taxon_summary = taxon_summary,
																	 verbose = verbose,
																	 criterion = criterion)
	message("DateLife search done!")
	return(res)
}

# print.datelife <- function(datelife){
# 	datelife_result_check(datelife$datelife_result)
# 	cat("Number and which queried taxa found in chronogram database")
# 	cat("Number of queried taxa not found")
# 	cat("Number of chronograms with at least two queried taxa found in database")
# }



#' Check if we obtained an empty search with the set of input taxon names
#' @inheritParams datelife_search
#' @param datelife_result An output of get_datelife_result function: A list of patristic matrices with names corresponding to original study citations.
#' @export
datelife_result_check <- function(datelife_result, use_tnrs = FALSE, verbose = FALSE){
	if(length(datelife_result) < 1) {
		warning("'datelife_result' is empty.", call. = FALSE)
		# if(verbose) {
			message("'input' taxon names were not found in the local chronogram database.")
			if(!use_tnrs) {
				message("setting 'use_tnrs = TRUE' might change this, but it is time consuming.")
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
				message("'datelife_result' has some elements that are not matrices; check it out.")
			}
		} else {
			message("'datelife_result' is not a list; check it out.")
		}
	}
	return(datelife_result)
}

#' \code{datelifeResult} object of three birds "Rhea americana", "Pterocnemia pennata", and "Struthio camelus"
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
