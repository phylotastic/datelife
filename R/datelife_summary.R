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
		median.result <- NULL
		overlap <- 2
		while(!inherits(median.result, "phylo")){
		  best_grove <- datelife::filter_for_grove(datelife_result,
		                criterion = "taxa", n = overlap)
		  median.result <- tryCatch(suppressMessages(datelife_result_median(best_grove)), error = function(e) NULL)
			# sometimes max(branching.times) is off (too big or too small), so we could
			# standardize by real median of original data (max(mrcas)).
			# median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
		  overlap <- overlap + 1
		}
		tree <- median.result
	}
	if(summary_format.in %in% c("newick_sdm", "phylo_sdm")) {
		sdm.result <- datelife_result_sdm(best_grove)
		tree <- sdm.result$phy
	}
	if(summary_format.in %in% c("newick_sdm", "newick_median")) {
		return.object <- ape::write.tree(tree)
	}
	if(summary_format.in %in% c("phylo_sdm", "phylo_median")) {
		return.object <- tree
	}
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
		# datelife_result <- best_grove
		unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
		good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose)
		if(length(good.matrix.indices) > 1) {
			unpadded.matrices <- unpadded.matrices[good.matrix.indices]
			SDM.result <- get_sdm(unpadded.matrices, weighting, verbose)
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
sdm_matrix_to_phylo <- function(sdm_matrix){
	# for testing:
	# datelife_result <- get_datelife_result(input = "cetacea")
    # unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
    # good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
    # if(length(good.matrix.indices) > 1) {
    #   unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    #   sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
    # }
	# which(sdm_matrix < 0)
	# sdm_matrix[ceiling(7301/ncol(sdm_matrix)),] # Eubalaena japonica,
	# sdm_matrix[,ceiling(261/nrow(sdm_matrix))]  # Eubalaena glacialis
	# xx <- sdm_matrix #[1:5, 1:5]
	#even removing negative values for small positive values gives back non ultrametric trees with njs
	# sdm_matrix[which(sdm_matrix < 0)] <- 0.01
	# test <- cluster_patristicmatrix(sdm_matrix)
	# class(test) <- "multiPhylo"
	# ape::is.ultrametric(test)
	# plot(test$njs)

	ages <- tA <- tB <- c() # compute the final length of the data frame: it's ncol(xx)^2 - sum(1:(ncol(xx)-1))
	# calibrations <- matrix(nrow = ncol(xx)^2 - sum(1:(ncol(xx)-1)), ncol = 3)
	# start <- ?
	for(i in seq(ncol(sdm_matrix))){
		# calibrations[start:start+i,1] <- xx[1:i,i]
		# calibrations[start:start+i,2] <- rownames(xx)[1:i]
		# calibrations[start:start+i,3] <- rep(colnames(xx)[i], i)
		# start <- sum(!is.na(calibrations[,1])) +1
		ages <- c(ages, sdm_matrix[1:i,i])
		tA <- c(tA, rownames(sdm_matrix)[1:i])
		tB <- c(tB, rep(colnames(sdm_matrix)[i], i))
	}
	calibrations <- data.frame(Age = ages, taxonA = tA, taxonB = tB)
	calibrations$taxonA <- as.character(calibrations$taxonA)
	calibrations$taxonB <- as.character(calibrations$taxonB)
	calibrations <- calibrations[!is.na(calibrations[,"Age"]), ] # get rid of NaN
	calibrations <- calibrations[calibrations[,"Age"] != 0, ] # get rid of 0's
	calibrations[calibrations[, "Age"] < 0, "Age"] <- 0.01 # replace negative values for a tiny number
	# chronogram <- geiger::PATHd8.phylo(phy_target, calibrations)
	# try(chronogram <- geiger::PATHd8.phylo(phy_target, calibrations), silent = TRUE)
	target_tree <- patristic_matrix_to_phylo(sdm_matrix)
	target_tree$edge.length <- NULL
	target_tree <- tree_add_nodelabels(tree = target_tree, node_index = "consecutive")  # all nodes need to be named so make_bladj_tree runs properly
	target_tree_nodes <- sapply(seq(nrow(calibrations)), function(i)
			phytools::findMRCA(tree = target_tree,
								 tips = as.character(calibrations[i,c("taxonA", "taxonB")]),
								 type = "node"))
	target_tree_nodes <- target_tree_nodes - ape::Ntip(target_tree)
	all_nodes <- sort(unique(target_tree_nodes))
	all_ages <- lapply(seq(length(all_nodes)), function(i) calibrations[target_tree_nodes == i, "Age"])
	# any(sapply(all_ages, is.null)) # all nodes have at least one calibration.
	calibrations2 <- data.frame(MRCA = paste0("n", all_nodes), MinAge = sapply(all_ages, min), MaxAge= sapply(all_ages, max), node = all_nodes)
	new.phy <- make_bladj_tree(tree = target_tree, nodenames = as.character(calibrations2$MRCA), nodeages = sapply(seq(nrow(calibrations2)), function(i) sum(calibrations2[i,c("MinAge", "MaxAge")])/2))
}
tree_add_dates2 <- function(calibrations, target_tree, method = "bladj"){
	if("bladj" %in% method){
		dated_tree_nodes <- sapply(seq(nrow(calibrations)), function(i)
				phytools::findMRCA(tree = target_tree,
									 tips = as.character(calibrations[i,c("taxonA", "taxonB")]),
									 type = "node"))
		dated_tree_nodes <- dated_tree_nodes - ape::Ntip(target_tree)
		# missing_taxa_phy$node.label[dated_tree_nodes] <- constraint_tree$calibrations$MRCA
		# cannot use hash number to name nodes, bladj collapses. So using "congNumber"
		target_tree$node.label[dated_tree_nodes] <- paste0("cong", seq(nrow(calibrations)))
		target_tree <- tree_add_nodelabels(tree = target_tree)  # all nodes need to be named so make_bladj_tree runs properly
		# this adds random names to unnamed nodes, but they have to coincide between target and reference
		# new.phy <- make_bladj_tree(tree = missing_taxa, nodenames = dated_tree$node.label, nodeages = tree_get_node_data(tree = dated_tree, node_data = "node_age")$node_age)
		new.phy <- make_bladj_tree(tree = target_tree, nodenames = target_tree$node.label[dated_tree_nodes], nodeages = sapply(seq(nrow(calibrations)), function(i) sum(calibrations[i,c("MinAge", "MaxAge")])/2))

	}
}
get_sdm <- function(unpadded.matrices, weighting, verbose){
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
	return(SDM.result)
}
get_goodmatrices <- function(unpadded.matrices, verbose){
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
	return(good.matrix.indices)
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
