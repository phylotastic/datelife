#' Get taxon summary of a datelifeResult object
#' @inheritParams datelife_query_check
#' @inheritParams datelife_result_check
#' @export
get_taxon_summary <- function(datelife_query = NULL, datelife_result){
	if(is.null(datelife_result) | !is.list(datelife_result)){
		return(NA)
		message("datelife_result argument must be a list of patristic matrices (you can get one with get_datelife_result()).")
	}
	if(is.null(datelife_query)){
		input <- NULL
		input.in <- unique(rapply(datelife_result, rownames))
		# if(taxon_summary.in != "none") {
			message("datelife_query argument is empty: showing taxon distribution of taxa found only in at least one chronogram. This excludes input taxa not found in any chronogram.")
		# }
	} else {
		# if(!is.character(input)) stop("input must be a character vector")
		input <- datelife_query_check(datelife_query = datelife_query)
		input.in <- input$cleaned_names
		# input.in <- input
	}
	# results.index <- datelife_result_study_index(datelife_result, cache)
	return.object <- NA
	input.match <- unique(rapply(datelife_result, rownames))
	# if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
	absent.input <- input.in[!input.in %in% input.match]
	if(length(absent.input) <= 0) {
		if(is.null(input)){
			absent.input <- NA # we cannot know if there are complete absent taxa because original query was not provided
		} else {
			absent.input <- "None" # we know there are no absent taxa
		}
	}
	taxon_list <- vector(mode = "list")
	# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
	for(result.index in sequence(length(datelife_result))){
		n <- rownames(datelife_result[[result.index]])
		m <- match(input.match,n)
		taxon_list[[result.index]] <- n[m]
	}
	taxon_matrix <- do.call(rbind, taxon_list) #transforms a list of names into a matrix of names
	taxon_matrix <- !is.na(taxon_matrix) # makes a boolean matrix
	colnames(taxon_matrix) <- input.match
	rownames(taxon_matrix) <- sequence(nrow(taxon_matrix))
	# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
	x <- rapply(datelife_result, rownames)
	prop <- c()
	for (taxon in input.match){
		prop <- c(prop, paste(length(which(taxon == x)), "/", length(datelife_result), sep=""))
	}
	taxon_summary <- data.frame(taxon = input.match, chronograms = prop)
	return(list(matrix = taxon_matrix, summary = taxon_summary, absent_taxa = absent.input))
}

#' Summarize a filtered results list from get_datelife_result function in various ways
#' @inheritParams datelife_query_check
#' @inheritParams datelife_result_check
#' @inheritParams datelife_search
#' @inherit datelife_search return details
#' @export
summarize_datelife_result <- function(datelife_query = NULL, datelife_result = NULL,
	summary_format = "phylo_all", partial = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"),
	summary_print = c("citations", "taxa"), taxon_summary = c("none", "summary", "matrix"),
	verbose = FALSE, criterion = c("trees", "taxa")) {
		# if(!partial) {
		# 	datelife_result <- datelife_result[which(!sapply(datelife_result, anyNA))]
		# } # not necessary cause already filtered in get_datelife_result
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	taxon_summ <- get_taxon_summary(datelife_query, datelife_result)
	if(length(taxon_summ) == 1){
		return(NA)
	}
	summary_format.in <- match.arg(summary_format, choices = c("citations", "mrca", "newick_all", "newick_sdm", "newick_median", "phylo_sdm", "phylo_median", "phylo_biggest", "phylo_all", "html", "data_frame"))
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
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- trees[which(!is.na(trees))]
		class(return.object) <- "multiPhylo"
	}
	if(summary_format.in == "phylo_biggest") {
		# enhance: choose from patristic matrix, not phylo objects, it will be faster.
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- get_biggest_phylo(trees) # NAs in trees are removed in get_biggest_phylo
	}
	# the following chunck is to test if n_overlap = 2 is enough to summarize results with sdm and median
	if(summary_format.in %in% c("newick_sdm", "phylo_sdm", "newick_median", "phylo_median")){
		best_grove <- get_best_grove(datelife_result, criterion = "taxa", n = 2)$best_grove
	}
	if(grepl("median", summary_format.in)){
		return.object <- datelife_result_median(best_grove)
	}
	if(grepl("sdm", summary_format.in)) {
		return.object <- datelife_result_sdm(best_grove)
	}
	if(grepl("newick", summary_format.in)) {
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
			if(summary_format.in !="data_frame") return.object <- c(return.object, list(taxon_distribution = taxon_summ$matrix))
		}
		if(taxon_summary.in == "summary") {
			return.object <- c(return.object, list(taxon_distribution = taxon_summ$summary))
		}
		return.object <- c(return.object, list(absent_taxa = data.frame(taxon = taxon_summ$absent_taxa)))
		names(return.object)[1] <- summary_format.in
		if(summary_format.in =="data_frame"){
			names(return.object)[1] <- "results"
			if(taxon_summary.in == "matrix") names(return.object)[1] <- "results_and_missing_taxa"
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

#' Get the tree with the most tips: the biggest tree
#' @param trees A list of trees as multiPhylo or as a plain list object.
#' @return A phylo object
#' @export
get_biggest_phylo <- function(trees){
	trees <- trees[which(!is.na(trees))] # removes NAs, which will return an error later on next logical:
	return.object <- trees[which(sapply(trees, ape::Ntip) == max(sapply(trees, ape::Ntip)))]
	if(length(return.object) >1 ) { #there are more than one tree with same number of taxa. Rather than take the first by default, take the one with the most intermediate depth (this assumes that the root node is the same for all trees). An example is the Bininda-Emonds et al mammal tree: there are three trees with min, max, and best guess calibrations. So, take the one in the middle.
		max.branching.time <- function(x) {
			return(max(ape::branching.times(x)))
		}
		tree.depths <- sapply(return.object, max.branching.time)
		return.object <- return.object[which.min(abs(tree.depths - stats::median(tree.depths)))]
	}
	if(class(return.object) != "phylo") {
		return.object <- return.object[[1]]
	}
	return.object
}
