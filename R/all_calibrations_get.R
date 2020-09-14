#' Get all calibrations from chronograms in a database (specified in cache).
#' @param input vector of names, a newick string, a phylo or multiPhylo object, a datelifeResult object
#' @param partial Boolean; default TRUE: use source trees even if they only match some of the desired taxa
#' @param use_tnrs Boolean; default False. If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match Boolean; default TRUE: use a slower TNRS to correct misspellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams datelife_search
#' @param each Boolean, default to FALSE (all calibrations are returned in the same data frame). If TRUE, calibrations from each chronogram are returned in a separate data frame.
#' @return A data_frame of calibrations
#' @export
# input <- tax_phyloallall[[1]]
get_all_calibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
		partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE,
		cache = get("opentree_chronograms"), verbose = FALSE, each = FALSE) {
	if(!inherits(input, "datelifeResult") & !inherits(input, "phylo") & !inherits(input, "multiPhylo")){
		datelife_phylo <- datelife_search(input = input, partial = partial, use_tnrs = use_tnrs,
			approximate_match = approximate_match, update_cache = update_cache, cache = cache,
			summary_format = "phylo_all", verbose = verbose)
	}
	# inherits(datelife_phylo, "datelifeResult")
	if(inherits(input, "datelifeResult")){
		datelife_phylo <- suppressMessages(summarize_datelife_result(datelife_result = input, summary_format = "phylo_all"))
	}
	if(inherits(input, "phylo")){
		datelife_phylo <- list(input)
	}
	if(inherits(input, "multiPhylo")){
		datelife_phylo <- input
	}
	if(each){
		constraints.df <- vector(mode= "list") # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
	} else {
		constraints.df <- data.frame() # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
	}
	for (i in seq(length(datelife_phylo))) {
		local.df <- suppressWarnings(geiger::congruify.phylo(reference = datelife_phylo[[i]],
			target = datelife_phylo[[i]], scale = NA, ncores = 1))$calibrations
		# suppressedWarnings bc of message when running geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA)
		# 		Warning message:
		# In if (class(stock) == "phylo") { :
		# the condition has length > 1 and only the first element will be used
		local.df$reference <- names(datelife_phylo)[i]
			if(each){
				constraints.df <- c(constraints.df, list(local.df))
			} else {
				if(i == 1) {
					constraints.df <- local.df
				} else {
					constraints.df <- rbind(constraints.df, local.df)
				}
		 }
	}
	if(each){
		names(constraints.df) <- names(datelife_phylo)
	}
	attr(constraints.df, "phylo_all") <- datelife_phylo
	return(constraints.df)
}
