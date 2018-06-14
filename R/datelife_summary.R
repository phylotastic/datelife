
#' Summarize a filtered results list from get_datelife_result function in various ways
#' @inheritParams datelife_query_check
#' @inheritParams datelife_result_check
#' @inheritParams datelife_search
#' @inherit datelife_search return details
#' @export
summarize_datelife_result <- function(datelife_query = NULL, datelife_result = NULL, summary_format = "phylo_all", partial = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), summary_print = c("citations", "taxa"), add_taxon_distribution = c("none", "summary", "matrix"), verbose = FALSE) {
		# if(!partial) {
		# 	datelife_result <- datelife_result[which(!sapply(datelife_result, anyNA))]
		# } # not necessary cause already filtered in get_datelife_result
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	if(is.null(datelife_result) | !is.list(datelife_result)){
		stop("datelife_result argument must be a list from get_datelife_result function.")
	}
	summary_format.in <- match.arg(summary_format, choices = c("citations", "mrca", "newick_all", "newick_sdm", "newick_median", "phylo_sdm", "phylo_median", "phylo_biggest", "phylo_all", "html", "data_frame"))
	add_taxon_distribution.in <- match.arg(add_taxon_distribution, choices = c("none", "summary", "matrix"))
	summary_print.in <- match.arg(summary_print, c("citations", "taxa", "none"), several.ok = TRUE)
	input <- datelife_query
	if(is.null(input)){
		input.in <- unique(rapply(datelife_result, rownames))
		# if(add_taxon_distribution.in != "none") {
			warning("datelife_query argument is empty: showing taxon distribution of taxa found only in at least one chronogram. This excludes input taxa not found in any chronogram.")
		# }
	} else {
		# if(!is.character(input)) stop("input must be a character vector")
		input <- datelife_query_check(datelife_query = input)
		input.in <- input$cleaned_names
		# input.in <- input
	}
	results.index <- datelife_result_study_index(datelife_result, cache)
	return.object <- NA
	input.match <- unique(rapply(datelife_result, rownames))
	# if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
	absent.input <- input.in[!input.in %in% input.match]
	if(length(absent.input) <= 0) {
		if(is.null(input)){
			absent.input <- "NULL"
		} else {
			absent.input <- "None"
		}
	}
	if(add_taxon_distribution.in == "matrix"){
		taxon_distribution_list <- vector(mode = "list")
		# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
		for(result.index in sequence(length(datelife_result))){
			n <- rownames(datelife_result[[result.index]])
			m <- match(input.match,n)
			taxon_distribution_list[[result.index]] <- n[m]
		}
		taxon_distribution_matrix <- do.call(rbind, taxon_distribution_list) #transforms a list of names into a matrix of names
		taxon_distribution_matrix <- !is.na(taxon_distribution_matrix) # makes a boolean matrix
		colnames(taxon_distribution_matrix) <- input.match
		rownames(taxon_distribution_matrix) <- sequence(nrow(taxon_distribution_matrix))
	}
	if(add_taxon_distribution.in == "summary" | any(grepl("taxa", summary_print.in))){ # may add here another condition: | makeup_brlen ==TRUE
		# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
		x <- rapply(datelife_result, rownames)
		prop <- c()
		for (taxon in input.match){
			prop <- c(prop, paste(length(which(taxon == x)), "/", length(datelife_result), sep=""))
		}
		taxon_distribution_summary <- data.frame(taxon = input.match, chronograms = prop)
	}
	if(summary_format.in == "citations") {
		return.object <- names(datelife_result)
	}
	if(summary_format.in == "mrca") {
		return.object <- datelife_result_MRCA(datelife_result, partial = partial)
	}
	if(summary_format.in == "newick_all") {
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "newick_sdm") {
		local.results <- datelife_result_sdm(datelife_result)
		datelife_result <- local.results$datelife_result
		tree <- local.results$phy
		return.object <- ape::write.tree(tree)
	}
	if(summary_format.in == "phylo_sdm") {
		local.results <- datelife_result_sdm(datelife_result)
		datelife_result <- local.results$datelife_result
		tree <- local.results$phy
		return.object <- tree
	}
	if(summary_format.in == "newick_median") {
		patristic.array <- patristic_matrix_list_to_array(datelife_result)
		median.matrix <- summary_patristic_matrix_array(patristic.array)
		tree <- patristic_matrix_to_newick(median.matrix)
		return.object <- tree
	}
	if(summary_format.in == "phylo_median") {
		patristic.array <- patristic_matrix_list_to_array(datelife_result)
		median.matrix <- summary_patristic_matrix_array(patristic.array)
		tree <- patristic_matrix_to_phylo(median.matrix)
		return.object <- tree
	}
	if(summary_format.in == "phylo_all") {
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "phylo_biggest") {
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- trees[which(sapply(trees, ape::Ntip)==max(sapply(trees, ape::Ntip)))]
		if(length(return.object)>1) { #there are more than one tree with same number of taxa. Rather than take the first by default, take the one with the most intermediate depth (this assumes that the root node is the same for all trees). An example is the Bininda-Emonds et al mammal tree: there are three trees with min, max, and best guess calibrations. So, take the one in the middle.
			max.branching.time <- function(x) {
				return(max(ape::branching.times(x)))
			}
			tree.depths <- sapply(return.object, max.branching.time)
			return.object <- return.object[which.min(abs(tree.depths - stats::median(tree.depths)))]
		}
		if(class(return.object)!="phylo") {
			return.object <- return.object[[1]]
		}
	}
	if(summary_format.in == "html") {
		out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
		if(add_taxon_distribution.in == "matrix"){
			out.vector1 <- paste(out.vector1, paste("</th><th>", colnames(taxon_distribution_matrix), sep="", collapse=""), sep="")
		}
		out.vector1 <- paste(out.vector1, "</th></tr>", sep="")
		ages <- datelife_result_MRCA(datelife_result, partial = partial)
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		out.vector2 <- c()
		for(result.index in sequence(length(datelife_result))) {
			out.vector2 <- paste(out.vector2, "<tr><td>",ages[result.index],"</td><td>",sum(!is.na(diag(datelife_result[[result.index]]))), "</td><td>", names(datelife_result)[result.index], "</td><td>", trees[result.index],  sep="")
			if(add_taxon_distribution.in == "matrix"){
				out.vector2 <- paste(out.vector2, paste("</td><td>", taxon_distribution_matrix[result.index,], sep="", collapse=""), sep="")
			}
			out.vector2 <- paste(out.vector2, "</td></tr>", sep="")
		}
		out.vector <- paste(out.vector1, out.vector2, "</table>", sep="")
		if(add_taxon_distribution.in == "summary"){
			taxon_distribution_summary.html <- as.matrix(taxon_distribution_summary)
			out.vector3 <- "<p></p><table border='1'><tr><th>taxon</th><th>chronograms</th><tr>"
			for (summary.index in sequence(nrow(taxon_distribution_summary.html))){
				out.vector3 <- paste(out.vector3, paste("</td><td>", taxon_distribution_summary.html[summary.index,], sep="", collapse=""), "</td></tr>", sep="")
			}
			out.vector <- paste(out.vector, out.vector3, "</table>", sep="")
		}
		# out.vector4 <- c()
		if(add_taxon_distribution.in != "none") {
			out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
			for (i in 1:length(absent.input)){
				out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", absent.input[i], "</td><tr>", sep="")
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
		if(add_taxon_distribution.in == "matrix"){
			out.df <- cbind(out.df, taxon_distribution_matrix)
		}
		rownames(out.df) <- NULL
		return.object <- out.df
	}
	if(add_taxon_distribution.in != "none" & summary_format.in !="html"){
		return.object <- list(return.object)
		if(add_taxon_distribution.in == "matrix") {
			if(summary_format.in !="data_frame") return.object <- c(return.object, list(taxon_distribution = taxon_distribution_matrix))
		}
		if(add_taxon_distribution.in == "summary") {
			return.object <- c(return.object, list(taxon_distribution = taxon_distribution_summary))
		}
		return.object <- c(return.object, list(absent_taxa = data.frame(taxon = absent.input)))
		names(return.object)[1] <- summary_format.in
		if(summary_format.in =="data_frame"){
			names(return.object)[1] <- "results"
			if(add_taxon_distribution.in == "matrix") names(return.object)[1] <- "results_and_missing_taxa"
		}
	}

	if(any("citations" %in% summary_print.in) & !any(summary_format.in %in% c("citations", "html", "data_frame"))) {
		if(summary_format.in == "citations"){
			message("Target taxa found in trees from:")
		} else {
			message("Source chronograms from:", "\n")
		}
		for (i in 1:length(datelife_result)){
			message(i, ": ", names(datelife_result)[i], "\n")
		}
	}
	if(any(grepl("taxa", summary_print.in)) & add_taxon_distribution.in!="summary") {
		message("Target taxa presence in source chronograms:")
		message(paste0(utils::capture.output(taxon_distribution_summary), collapse = "\n"), "\n")
		message("Target taxa completely absent from source chronograms:")
		message(paste0(utils::capture.output(data.frame(taxon = absent.input)), collapse = "\n"), "\n")
	}
	return(return.object)
}

#' Function to compute the SDM supertree (Criscuolo et al. 2006) from a datelifeResult object.
#' @inheritParams datelife_result_check
#' @param weighting A character vector indicating how much weight to give to each input tree in the SDM analysis.
#' 	 Choose one of:
#' \describe{
#'	\item{weighting = "flat"}{All trees have equal weighting.
#'	}
#'	\item{weighting = "taxa"}{Weight is proportional to number of taxa.
#'	}
#'	\item{weighting = "inverse"}{Weight is proportional to 1 / number of taxa.
#'	}
#' }
#' @return A list of two elements:
#' \describe{
#'	\item{phy}{A phylo object with the output chronogram from SDM analysis.
#'	}
#'	\item{data}{A datelifeResult object with the chronograms that were used to construct the SDM tree.
#'	}
#' }
#' @export
#' @details
#' Criscuolo A, Berry V, Douzery EJ, Gascuel O. SDM: a fast distance-based approach for (super) tree building in phylogenomics. Syst Biol. 2006. 55(5):740. doi: 10.1080/10635150600969872.
datelife_result_sdm <- function(datelife_result, weighting = "flat", verbose = TRUE) {
	phy <- NA
	used.studies <- names(datelife_result)
	unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
	good.matrix.indices <- c()
	for(i in sequence(length(unpadded.matrices))) {
		test.result <- NA
		# Rationale here: some chronograms always cause errors with SDM, even when trying to get a consensus of them
		# with themselves. For now, throw out of synthesis.
		try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent = TRUE)
		if (verbose){
			message(cat(i, "out of", length(unpadded.matrices), "chronograms tried: "), appendLF = FALSE)
		}
		if(is.finite(test.result)) {
			good.matrix.indices <- append(good.matrix.indices,i)
			if (verbose){
				message(cat(" Ok."))
			}
		} else {
			if (verbose){
				message(cat(" Failed."))
			}
		}
	}
	if(length(good.matrix.indices) > 0) {
		unpadded.matrices <- unpadded.matrices[good.matrix.indices]
		used.studies <- used.studies[good.matrix.indices]
		weights = rep(1, length(unpadded.matrices))
		if (weighting=="taxa") {
			weights = unname(sapply(unpadded.matrices, dim)[1,])
		}
		if (weighting=="inverse") {
			weights = 1/unname(sapply(unpadded.matrices, dim)[1,])
		}
		if (verbose){
			message(cat("\n", "Synthesizing", length(unpadded.matrices), "chronograms with SDM"))
		}
		SDM.result <- do.call(ape::SDM, c(unpadded.matrices, weights))[[1]]
		phy <- patristic_matrix_to_phylo(SDM.result, clustering_method = "upgma")  # nj or njs do not work if patristic matrices output from sdm. Go back to this later

	} else {
		warning("All input chronograms throw an error when running SDM. This is not your fault.")
		stop("SDM cannot be run with this set of chronograms.")
	}
	class(unpadded.matrices) <- "datelifeResult"
	return(list(phy = phy, data = unpadded.matrices))
}

#' Get an otol induced dated subtree from your set of queried taxa
#' @inheritParams datelife_search
#' @export
#' @details otol dated tree from Stephen Smith's http://141.211.236.35:10999/
get_dated_otol_induced_subtree <- function(input){
	input <- datelife_query_check(input)
	input_ott_match <- rotl::tnrs_match_names(input$cleaned_names)$ott_id
  	system(paste0('curl -X POST http://141.211.236.35:10999/induced_subtree -H "content-type:application/json" -d \'{"ott_ids":[', paste(input_ott_match, collapse = ", "), "]}'", "> /tmp/otol_induced_subtree.txt"))
  	otol_induced_subtree <- jsonlite::fromJSON(txt = "/tmp/otol_induced_subtree.txt")
  	system(paste0('curl -X POST http://141.211.236.35:10999/rename_tree -H "content-type:application/json" -d \'{"newick":"', otol_induced_subtree$newick, "\"}'", "> /tmp/otol_induced_subtree_names.txt"))
  	otol_induced_subtree_names <- jsonlite::fromJSON(txt = "/tmp/otol_induced_subtree_names.txt")
  	phy <- ape::read.tree(text = otol_induced_subtree_names$newick)
  	return(phy)
}
