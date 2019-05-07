# get_all_calibrations.default method applies to character vector
# it gets the dateliferesult object (and discriminates between vector of names and newick string)
# then calls get_all_calibrations.dateliferesult
# which calls get_all_calibrations.phylo or get_all_calibrations.multiphylo or
# this one gets the actual calibrations with congruify
# the other method is get_all_calibrations.data.frame alias get_all_calibrations.matrix
# that uses the code in summarize_summary_matrix

# #' Get all calibrations from a vector of names, a tree as newick string or phylo, a list of trees, or a patristic distance matrix.
# #' @param input A a phylo or multiphylo object, or a patristic distance matrix with taxon names as column and row names.
# #' @return An data frame with calibrations.
# #' @export
# get_all_calibrations <- function(input) {
#     UseMethod("get_all_calibrations", input)
# }

# #' @return \code{NULL}
# #'
# #' @rdname get_all_calibrations
# #' @method get_all_calibrations default

# get_all_calibrations.multiphylo <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),

#' Get all calibrations from chronograms in a database (specified in cache).
#' @param input vector of names, a newick string, a phylo or multiPhylo object, a datelifeResult object
#' @param partial Boolean; default TRUE: use source trees even if they only match some of the desired taxa
#' @param use_tnrs Boolean; default False. If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match Boolean; default TRUE: use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams datelife_search
#' @return A data_frame of calibrations
#' @export
get_all_calibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
		partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE,
		cache = get("opentree_chronograms"), verbose = FALSE) {
	if(!inherits(input, "datelifeResult") & !inherits(input, "phylo") & !inherits(input, "multiPhylo")){
		datelife_phylo <- datelife_search(input = input, partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = update_cache, cache = cache, summary_format = "phylo_all", verbose = verbose)
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
	constraints.df <- data.frame() # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
	for (i in seq(length(datelife_phylo))) {
		local.df <- suppressWarnings(geiger::congruify.phylo(reference = datelife_phylo[[i]],
			target = datelife_phylo[[i]], scale = NA, ncores = 1))$calibrations
		# suppressedWarnings bc of meesage when running geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA)
		# 		Warning message:
		# In if (class(stock) == "phylo") { :
		# the condition has length > 1 and only the first element will be used
		local.df$reference <- names(datelife_phylo)[i]
		if(i == 1) {
			constraints.df <- local.df
		} else {
			constraints.df <- rbind(constraints.df, local.df)
		}
	}
	return(constraints.df)
}

#' Internal for summary_matrix_to_phylo().
#' @inheritParams summary_matrix_to_phylo
summarize_summary_matrix <- function(summ_matrix){
      	ages <- tA <- tB <- c()
        # to compute the final length of the data frame do: ncol(xx)^2 - sum(1:(ncol(xx)-1))
    	# calibrations <- matrix(nrow = ncol(xx)^2 - sum(1:(ncol(xx)-1)), ncol = 3)
    	# identify if SDM matrix has some negative values; extract taxon names:
    	negs <- which(summ_matrix < 0)
    	neg_names <- rownames(summ_matrix)[ceiling(negs/nrow(summ_matrix))]
    	# extract unique ages from summ_matrix:
    	for(i in seq(ncol(summ_matrix))){
    		ages <- c(ages, summ_matrix[1:i,i])
    		tA <- c(tA, rownames(summ_matrix)[1:i])
    		tB <- c(tB, rep(colnames(summ_matrix)[i], i))
    	}
        # tA <- gsub(" ", "_", tA)
        # tB <- gsub(" ", "_", tB)
    	calibrations <- data.frame(Age = ages, taxonA = tA, taxonB = tB, stringsAsFactors = FALSE)
    	calibrations <- calibrations[!is.na(calibrations[,"Age"]), ] # get rid of NaN and NAs
    	calibrations <- calibrations[calibrations[,"Age"] != 0, ] # get rid of 0's
    	calibrations <- calibrations[calibrations[,"Age"] > 0, ] # get rid of negative values too
    	if(any(is.na(calibrations[,"Age"]))){
    		warning("for some reason there are still NAs in the matrix")}
    	# enhance: where does this negative values come from in SDM?
        return(calibrations)
}
