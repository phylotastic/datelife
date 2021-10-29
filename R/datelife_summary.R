#' Get a taxon summary of a \code{datelifeResult} object
#'
# #' To be renamed summary_taxon.
#'
#' @param datelife_query A \code{datelifeQuery} object, output of \code{\link[=make_datelife_query]{make_datelife_query}}.
#' @inheritParams summarize_datelife_result
#' @export
get_taxon_summary <- function(datelife_result = NULL,
															datelife_query = NULL){

	# datelife_result <- check_datelife_result(datelife_result)
	if (is.null(datelife_result) | !inherits(datelife_result, "datelifeResult")) {
		warning("'datelife_result' argument must be a list of patristic matrices (you can get one with get_datelife_result()).")
		return(NA)
	}

	if (suppressMessages(is_datelife_query(datelife_query))) {
		if (is.null(attributes(datelife_result)$datelife_query)) {
			cleaned_names <- datelife_query$cleaned_names
		} else {
			input <- attributes(datelife_result)$datelife_query
			cleaned_names <- attributes(datelife_result)$datelife_query$cleaned_names
		}
	} else {
		message("'datelife_query' argument was not provided.")
		message("Taxa absent in all chronograms are not reported.")
		cleaned_names <- unique(rapply(datelife_result, rownames))
	}
	# results.index <- datelife_result_study_index(datelife_result, cache)
	return.object <- NA
	input.match <- unique(rapply(datelife_result, rownames))
	# if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
	absent.input <- cleaned_names[!cleaned_names %in% input.match]
	if (length(absent.input) <= 0) {
		if (is.null(datelife_query)) {
			absent.input <- NA # we cannot know if there are any taxa that are completely absent, because original query was not provided
		} else {
			absent.input <- "None" # we know there are no absent taxa
		}
	}
	taxon_list <- vector(mode = "list")
	# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
	for (result.index in sequence(length(datelife_result))) {
		n <- rownames(datelife_result[[result.index]])
		m <- match(input.match,n)
		taxon_list[[result.index]] <- n[m]
	}
	taxon_matrix <- do.call(rbind, taxon_list) #transforms a list of names into a matrix of names
	taxon_matrix <- !is.na(taxon_matrix) # makes a boolean matrix
	colnames(taxon_matrix) <- input.match
	rownames(taxon_matrix) <- paste0("Chronogram", sequence(nrow(taxon_matrix)))
	chronogram_names <- names(datelife_result)
	names(chronogram_names) <- rownames(taxon_matrix)
	# from here, replace by print: take it to a print function, so we only store the matrix, absent taxa and chronogram names
	# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
	x <- rapply(datelife_result, rownames)
	prop <- c()
	for (taxon in input.match) {
		prop <- c(prop, paste0(length(which(taxon == x)), "/", length(datelife_result)))
	}
	taxon_summary <- data.frame(taxon = input.match, chronograms = prop)
	taxon_number <- sapply(seq(nrow(taxon_matrix)), function(x) sum(taxon_matrix[x,]))
	taxon_summary2 <- data.frame(chronogram = names(datelife_result),
		taxon_number = taxon_number, total_taxa = rep(ncol(taxon_matrix), length(datelife_result)))
	res <- list(matrix = taxon_matrix, summary = taxon_summary, summary2 = taxon_summary2,
		absent_taxa = unique(absent.input))
	# end replace by print
	# res <- list(matrix = taxon_matrix, chronogram_names = chronogram_names, absent_taxa = unique(absent.input))
	class(res) <- "datelifeTaxonSummary"
	return(res)
}
# print.datelifeTaxonSummary <- function(taxon_summary){
#
# }
#' Summarize a \code{datelifeResult} object
#'
#' @description Get different types of summaries from a \code{datelifeResult}
#' object, an output from [get_datelife_result()].
#' This allows rapid processing of data.
#' If you need a list of chronograms from your \code{datelifeResult} object, this
#' is the function you are looking for.
#'
#' @param datelife_result A \code{datelifeResult} object, an output
#' of [get_datelife_result()].
#' @inheritParams get_taxon_summary
#' @inheritParams datelife_search
#' @inherit datelife_search return details
#' @export
summarize_datelife_result <- function(datelife_result = NULL,
																			datelife_query = NULL,
																			summary_format = "phylo_all",
																			partial = TRUE,
																			update_opentree_chronograms = FALSE,
																			cache = "opentree_chronograms",
																			summary_print = c("citations", "taxa"),
																			taxon_summary = c("none", "summary", "matrix"),
																			verbose = FALSE,
																			criterion = c("trees", "taxa")) {
	if(update_opentree_chronograms){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	} else {
		if("opentree_chronograms" %in% cache){
			utils::data("opentree_chronograms", package = "datelife")
			cache <- get("opentree_chronograms")
		}
	}
	taxon_summ <- get_taxon_summary(datelife_result = datelife_result,
																	datelife_query = datelife_query)
	if(length(taxon_summ) == 1){
		message("get_taxon_summary failed.")
		return(NA)
	}
	summary_format.in <- match.arg(summary_format, choices = c("citations",
																														 "mrca",
																														 "newick_all",
																														 "newick_sdm",
																														 "newick_median",
																														 "phylo_sdm",
																														 "phylo_median",
																														 "phylo_biggest",
																														 "phylo_all",
																														 "html",
																														 "data_frame"))
	taxon_summary.in <- match.arg(taxon_summary, choices = c("none", "summary", "matrix"))
	summary_print.in <- match.arg(summary_print, c("citations", "taxa", "none"), several.ok = TRUE)
	if(summary_format.in == "citations") {
		return.object <- names(datelife_result)
	}
	mrcas <- datelife_result_MRCA(datelife_result, partial = partial)  # this is later used for median and sdm
	if(summary_format.in == "mrca") {
		return.object <- mrcas
	}
	if(summary_format.in == "newick_all") {
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "phylo_all") {
		trees <- suppressWarnings(lapply(datelife_result, patristic_matrix_to_phylo)) # suppress warning "Converting from patristic distance matrix to a tree resulted in some negative branch lengths"
		return.object <- trees[which(!is.na(trees))]
		class(return.object) <- "multiPhylo"
	}
	if(summary_format.in == "phylo_biggest") {
		# enhance: choose from patristic matrix, not phylo objects, it will be faster.
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- get_biggest_phylo(trees) # NAs in trees are removed in get_biggest_phylo
	}
	# test if n_overlap = 2 is enough to summarize results with sdm and median:
	if(summary_format.in %in% c("newick_sdm", "phylo_sdm", "newick_median", "phylo_median")){
		best_grove <- get_best_grove(datelife_result, criterion = "taxa", n = 2)$best_grove
	}
	if(inherits(datelife_query, "datelifeQuery")){
		target_tree <- datelife_query$phy
	} else {
		target_tree <- NULL
	}
	if(grepl("median", summary_format.in)){
		return.object <- datelife_result_median(best_grove, target_tree = target_tree)
	}
	if(grepl("sdm", summary_format.in)) {
		return.object <- datelife_result_sdm(best_grove, target_tree = target_tree)
	}
	if(summary_format.in %in% c("newick_sdm", "newick_median")) {
		return.object <- ape::write.tree(return.object)
	}
	# if(summary_format.in %in% c("phylo_sdm", "phylo_median")) {
	# 	return.object <- tree
	# }
	if(summary_format.in == "html") {
		out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
		if(taxon_summary.in == "matrix"){
			out.vector1 <- paste(out.vector1, paste("</th><th>", colnames(taxon_summ$matrix), sep="", collapse=""), sep="")
		}
		out.vector1 <- paste(out.vector1, "</th></tr>", sep="")
		ages <- datelife_result_MRCA(datelife_result, partial = partial)
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		out.vector2 <- c()
		for(result.index in sequence(length(datelife_result))) {
			out.vector2 <- paste(out.vector2, "<tr><td>",ages[result.index],"</td><td>",sum(!is.na(diag(datelife_result[[result.index]]))), "</td><td>", names(datelife_result)[result.index], "</td><td>", trees[result.index],  sep="")
			if(taxon_summary.in == "matrix"){
				out.vector2 <- paste(out.vector2, paste("</td><td>", taxon_summ$matrix[result.index,], sep="", collapse=""), sep="")
			}
			out.vector2 <- paste(out.vector2, "</td></tr>", sep="")
		}
		out.vector <- paste(out.vector1, out.vector2, "</table>", sep="")
		if(taxon_summary.in == "summary"){
			taxon_summ$summary.html <- as.matrix(taxon_summ$summary)
			out.vector3 <- "<p></p><table border='1'><tr><th>taxon</th><th>chronograms</th><tr>"
			for (summary.index in sequence(nrow(taxon_summ$summary.html))){
				out.vector3 <- paste(out.vector3, paste("</td><td>", taxon_summ$summary.html[summary.index,], sep="", collapse=""), "</td></tr>", sep="")
			}
			out.vector <- paste(out.vector, out.vector3, "</table>", sep="")
		}
		if(taxon_summary.in != "none") {
			out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
			for (i in 1:length(taxon_summ$absent_taxa)){
				out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", taxon_summ$absent_taxa[i], "</td><tr>", sep="")
			}
			out.vector4 <- paste(out.vector4, "</table>", sep="")
			out.vector <- paste(out.vector, out.vector4, sep="")
		}
		return.object <- out.vector
	}
	if(summary_format.in == "data_frame") {
		out.df <- data.frame()
		ages <- datelife_result_MRCA(datelife_result, partial = partial)
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		for(result.index in sequence(length(datelife_result))) {
			out.line<- data.frame(Age = ages[result.index],Ntax = sum(!is.na(diag(datelife_result[[result.index]]))), Citation = names(datelife_result)[result.index], Newick= trees[result.index])
			if(result.index == 1) {
				out.df <- out.line
			} else {
				out.df <- rbind(out.df, out.line)
			}
		}
		if(taxon_summary.in == "matrix"){
			out.df <- cbind(out.df, taxon_summ$matrix)
		}
		rownames(out.df) <- NULL
		return.object <- out.df
	}
	if(taxon_summary.in != "none" & summary_format.in !="html"){
		return.object <- list(return.object)
		if(taxon_summary.in == "matrix") {
			if(summary_format.in !="data_frame"){
				return.object <- c(return.object, list(taxon_distribution = taxon_summ$matrix))}
		}
		if(taxon_summary.in == "summary") {
			return.object <- c(return.object, list(taxon_distribution = taxon_summ$summary))
		}
		return.object <- c(return.object, list(absent_taxa = data.frame(taxon = taxon_summ$absent_taxa)))
		names(return.object)[1] <- summary_format.in
		if(summary_format.in =="data_frame"){
			names(return.object)[1] <- "results"
			if(taxon_summary.in == "matrix") {
				names(return.object)[1] <- "results_and_missing_taxa"}
		}
	}
	if(any("citations" %in% summary_print.in) & !any(summary_format.in %in% c("citations", "html", "data_frame"))) {
		if(summary_format.in == "citations"){
			message("Input taxa found in trees from:")
		} else {
			message("Source chronograms from:", "\n")
		}
		for (i in 1:length(datelife_result)){
			message(i, ": ", names(datelife_result)[i], "\n")
		}
	}
	if(any(grepl("taxa", summary_print.in)) & taxon_summary.in!="summary") {
		message("Input taxa presence across source chronograms:")
		message(paste0(utils::capture.output(taxon_summ$summary), collapse = "\n"), "\n")
		message("Input taxa completely absent from source chronograms:")
		message(paste0(utils::capture.output(data.frame(taxon = taxon_summ$absent_taxa)), collapse = "\n"), "\n")
	}
	return(return.object)
}
# TODO: method to summarize a dateliferesult object
# #' Main function to summarize a datelifeResult object
# #' @param object A "datelifeResult" object, typically an output of \code{\link[=get_datelife_result]{get_datelife_result}.
# #' @param ... further arguments passed to or from other methods
# #' @param partial Boolean for whether to include partial matches
# #' @method summary datelifeResult
# #' @export
# summary.datelifeResult <- function(object, ..., partial = TRUE){
# 	mrcas <- datelife_result_MRCA(object, partial = partial)
# 	res <- list(mrca = mrcas)
# 	class(res) <- "datelifeResultSummary"
# 	return(res)
# }

#' Get the tree with the most tips: the biggest tree
#'
#'
#'
#' @param trees A list of trees as multiPhylo or as a plain list object.
#' @return A phylo object with a citation slot with the citation of the biggest tree
#' @export
get_biggest_phylo <- function(trees){
	trees <- trees[which(!is.na(trees))] # removes NAs, which will return an error later on next logical:
	tree_citation <- names(trees)
	return.object <- trees[which(sapply(trees, ape::Ntip) == max(sapply(trees, ape::Ntip)))]
	tree_citation <- tree_citation[which(sapply(trees, ape::Ntip) == max(sapply(trees, ape::Ntip)))]
	if(length(return.object) >1 ) {  # there are more than one tree with same number of taxa. Rather than take the first by default, take the one with the most intermediate depth (this assumes that the root node is the same for all trees). An example is the Bininda-Emonds et al. mammal tree: there are three trees with min, max, and best guess calibrations. So, take the one in the middle.
		max.branching.time <- function(x) {
			return(max(ape::branching.times(x)))
		}
		tree.depths <- sapply(return.object, max.branching.time)
		return.object <- return.object[which.min(abs(tree.depths - stats::median(tree.depths)))]
		tree_citation <- tree_citation[which.min(abs(tree.depths - stats::median(tree.depths)))]
	}
	if(class(return.object) != "phylo") {
		return.object <- return.object[[1]]
	}
	return.object$citation <- tree_citation
	return.object
}
