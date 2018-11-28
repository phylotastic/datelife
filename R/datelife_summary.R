
#' Summarize a filtered results list from get_datelife_result function in various ways
#' @inheritParams datelife_query_check
#' @inheritParams datelife_result_check
#' @inheritParams datelife_search
#' @inherit datelife_search return details
#' @export
summarize_datelife_result <- function(datelife_query = NULL, datelife_result = NULL,
	summary_format = "phylo_all", partial = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"),
	summary_print = c("citations", "taxa"), add_taxon_distribution = c("none", "summary", "matrix"),
	verbose = FALSE, criterion = c("trees", "taxa")) {
		# if(!partial) {
		# 	datelife_result <- datelife_result[which(!sapply(datelife_result, anyNA))]
		# } # not necessary cause already filtered in get_datelife_result
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	if(is.null(datelife_result) | !is.list(datelife_result)){
		stop("datelife_result argument must be a list of patristic matrices (you can get one with get_datelife_result()).")
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
	mrcas <- datelife_result_MRCA(datelife_result, partial = partial)  # this is later used for median and sdm
	if(summary_format.in == "mrca") {
		return.object <- mrcas
	}
	if(summary_format.in == "newick_all") {
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "phylo_all") {
		trees <- lapply(datelife_result, patristic_matrix_to_phylo, clustering_method = "nj")
		return.object <- trees[which(!is.na(trees))]
		class(return.object) <- "multiPhylo"
	}
	if(summary_format.in == "phylo_biggest") {
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- get_biggest_phylo(trees)
	}
	# the following chunck is to test if n_overlap = 2 is enough to summarize results with sdm and median
	if(summary_format.in %in% c("newick_sdm", "phylo_sdm", "newick_median", "phylo_median")){
		median.result <- NULL
		overlap <- 2
		while(is.null(median.result)){
		  best_grove <- datelife::filter_for_grove(datelife_result,
		                criterion = "taxa", n = overlap)
		  median.result <- tryCatch(datelife_result_median(best_grove), error = function(e) NULL)
			# sometimes max(branching.times) is off (too big or too small), so we could
			# standardize by real median of original data (max(mrcas)).
			# median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
		  overlap <- overlap + 1
		}
		tree <- median.result
	}
	if(summary_format.in %in% c("newick_sdm", "phylo_sdm")) {
		#sometimes the best_grove for median does not work for sdm, so we tryCatch te result:
		# sdm.result <- tryCatch(datelife_result_sdm(best_grove), error = function(e) NULL)
		# while(is.null(sdm.result)){
		# 	best_grove <- datelife::filter_for_grove(datelife_result,
		#                 criterion = "taxa", n = overlap)  # we try the last overlap value tried from median
		#   sdm.result <- tryCatch(datelife_result_sdm(best_grove), error = function(e) NULL)
		# 	# sometimes max(branching.times) is off (too big or too nsmall), so  we
		# 	# standardize by real median of original data (max(mrcas)).
		# 	# median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
		#   overlap <- overlap + 1
		# }
		# datelife_result <- sdm.results$data
		sdm.result <- datelife_result_sdm(best_grove)
		tree <- sdm.result$phy
		# rm(sdm.result, best_grove)
	}
	if(summary_format.in %in% c("newick_sdm", "newick_median")) {
		return.object <- ape::write.tree(tree)
	}
	if(summary_format.in %in% c("phylo_sdm", "phylo_median")) {
		return.object <- tree
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

#' Function to compute median of a datelifeResult object.
#' @inheritParams patristic_matrix_to_phylo
#' @inheritParams datelife_result_check
datelife_result_median <- function(datelife_result, clustering_method = "nj") {
	patristic.array <- patristic_matrix_list_to_array(datelife_result)
	median.matrix <- summary_patristic_matrix_array(patristic.array)
	# when matrix comes from median, upgma gives much older ages than expected
	# we use nj to cluster in this case
	median.phylo <- patristic_matrix_to_phylo(median.matrix, clustering_method = clustering_method, fix_negative_brlen = TRUE)
	return(median.phylo)
}

#' Function to compute the SDM supertree (Criscuolo et al. 2006) from a datelifeResult object.
#' @inheritParams datelife_result_check
#' @inheritParams patristic_matrix_to_phylo
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
datelife_result_sdm <- function(datelife_result, weighting = "flat", verbose = TRUE, clustering_method = "nj") {
	# add check datelife_result
	phy <- NA
	# used.studies <- names(datelife_result)
	if(length(datelife_result) == 1){
		phy <- patristic_matrix_to_phylo(datelife_result, clustering_method = clustering_method, fix_negative_brlen = FALSE)
		unpadded.matrices <- datelife_result
	} else {
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
		if(length(good.matrix.indices) > 1) {
			unpadded.matrices <- unpadded.matrices[good.matrix.indices]
			# used.studies <- used.studies[good.matrix.indices]
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
			# it is important to use upgma as clustering method; nj produces much younger ages when the matrix comes from sdm
			phy <- patristic_matrix_to_phylo(SDM.result, clustering_method = clustering_method, fix_negative_brlen = TRUE)

		} else {
			if(length(good.matrix.indices) == length(datelife_result)) {
				warning("There are not enough input chronograms to run SDM. You can help uploading trees to Open Tree of Life tree store.")
			} else {
				warning("All input chronograms throw an error when running SDM. This is not your fault.")
			}
			stop("SDM tree cannot be generated with this set of taxa.")
		}
		class(unpadded.matrices) <- "datelifeResult"
	}
	return(list(phy = phy, data = unpadded.matrices))
}

#' Get an otol induced dated subtree from your set of queried taxa
#' @inheritParams datelife_search
#' @param ott_ids Numeric vector of Open Tree Taxonomy ids
#' @return A phylo object with edge length proportional to time in Myrs. NA if 1 or no inputs are valid.
#' @export
#' @details otol dated tree from Stephen Smith's otol scaling service
get_dated_otol_induced_subtree <- function(input = c("Felis silvestris", "Homo sapiens"), ott_ids = NULL){
	if(is.null(ott_ids)){
		input <- datelife_query_check(input)
		input_ott_match <- suppressWarnings(as.numeric(rotl::tnrs_match_names(input$cleaned_names)$ott_id))
		if(any(is.na(input_ott_match))){
			message(paste0("Input '", paste(input[which(is.na(input_ott_match))], collapse = "', '"), "', not found in Open Tree of Life Taxonomy."))
			input_ott_match <- input_ott_match[which(!is.na(input_ott_match))]
		}
	} else {
		input_ott_match <- suppressWarnings(as.numeric(ott_ids))
		if(any(is.na(input_ott_match))){
			message(paste0("Ott ids '", paste(ott_ids[which(is.na(input_ott_match))], collapse = "', '"), "', not numeric and will be excluded from the search."))
			input_ott_match <- input_ott_match[which(!is.na(input_ott_match))]
		}
	}
	if(length(input_ott_match) < 2){
		message("At least two valid names or numeric ott_ids are needed to get a tree")
		return(NA)
	}
  pp <- tryCatch(httr::POST("http://141.211.236.35:10999/induced_subtree",
						body = list(ott_ids = input_ott_match),
						encode = "json", httr::timeout(10)), error = function(e) NA)
	if(length(pp) > 1){
		pp <- httr::content(pp)
		rr <- httr::POST("http://141.211.236.35:10999/rename_tree",
		          body = list(newick = pp$newick),
		          encode = "json", httr::timeout(10))
		rr <- httr::content(rr)
		rr <- ape::read.tree(text = rr$newick)
		rr$ott_ids <- ape::read.tree(text = pp$newick)$tip.label
	  return(rr)
	}	else {
		return(pp)
	}
}

#' Get the tree with the most tips: the biggest tree
#' @param trees A list of trees as multiPhylo or as a plain list object.
#' @return A phylo object
#' @export
get_biggest_phylo <- function(trees){
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
