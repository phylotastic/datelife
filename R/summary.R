#' Summarize a `datelifeResult` object.
#' @param object An object of class `datelifeResult`, usually an output of [get_datelife_result()].
#' @inheritParams get_taxon_summary
#' @param na_rm Default to `TRUE`, whether to include partial matches or not.
#' @param ... Further arguments passed to or from other methods.
#' @return A named `list` of 11 elements:
#' \describe{
#' 	 \item{"citations"}{A character vector of references where chronograms with
#' 	 				some or all of the target taxa are published (source chronograms).}
#' 	 \item{"mrca"}{A named numeric vector of most recent common ancestor (mrca)
#' 	 				ages of target taxa defined in input, obtained from the source chronograms.
#' 	 				Names of mrca vector are equal to citations.}
#' 	 \item{"newick_all"}{A named character vector of newick strings corresponding
#' 	 				to target chronograms derived from source chronograms. Names of newick_all
#' 	 				vector are equal to citations.}
#' 	 \item{"newick_sdm"}{Only if multiple source chronograms are available. A
#' 	 				character vector with a single newick string corresponding to a target
#' 	 				chronogram obtained with SDM supertree method (Criscuolo et al. 2006).}
#' 	 \item{"newick_median"}{Only if multiple source chronograms are available.
#' 	 				A character vector with a single newick string corresponding to a target
#' 	 				chronogram from the median of all source chronograms.}
#' 	 \item{"phylo_sdm"}{Only if multiple source chronograms are available. A
#' 	 				phylo object with a single target chronogram obtained with SDM supertree
#' 	 				method (Criscuolo et al. 2006).}
#' 	 \item{"phylo_median"}{Only if multiple source chronograms are available. A
#' 	 				phylo object with a single target chronogram obtained from source
#' 	 				chronograms with median method.}
#' 	 \item{"phylo_all"}{A named list of phylo objects corresponding to each target
#' 	 				chronogram obtained from available source chronograms. Names of
#' 	 				phylo_all list correspond to citations.}
#' 	 \item{"phylo_biggest"}{The chronogram with the most taxa. In the case of a
#' 	 				tie, the chronogram with clade age closest to the median age of the
#' 	 				equally large trees is returned.}
#' 	 \item{"html"}{A character vector with an html string that can be saved and
#' 	 				then opened in any web browser. It contains a 4 column table with data on
#' 	 				target taxa: mrca, number of taxa, citations of source chronogram and
#' 	 				newick target chronogram.}
#' 	 \item{"data_frame"}{A 4 column `data.frame` with data on target taxa: mrca, number of
#' 	 				taxa, citations of source chronograms and newick string.}
#' }
#' @export
summary.datelifeResult <- function(object,
	                                 datelife_query,
																	 na_rm = TRUE,
																   ...) {
	datelife_result <- object
	mrcas <- datelife_result_MRCA(datelife_result, na_rm = na_rm)
	newick_all <- datelife_result_newick_all(datelife_result, na_rm = na_rm)
	phylo_all <- datelife_result_phylo_all(datelife_result, na_rm = na_rm)
	biggest <- get_biggest_multiphylo(phylo_all) # NAs in trees are removed in get_biggest_phylo

	median_and_sdm <- datelife_result_median_and_sdm(datelife_result, datelife_query, na_rm = na_rm)
  html_table <- datelife_result_html(datelife_result, datelife_query, na_rm = na_rm)
  data_frame <- datelife_result_data_frame(datelife_result, na_rm = na_rm)
	res <- list(citations = names(mrcas),
	            mrca = mrcas,
		          newick_all = newick_all,
							phylo_all = phylo_all,
						  newick_median = median_and_sdm$median_newick,
						  phylo_median = median_and_sdm$median_phylo,
							newick_sdm = median_and_sdm$sdm_newick,
							phylo_sdm = median_and_sdm$sdm_phylo,
							phylo_biggest = biggest,
							html = html_table,
							data_frame = data_frame)
	class(res) <- "datelifeResultSummary"
	return(res)
}

#' Get a numeric vector of MRCAs from a `datelifeResult` object. Used in [summarize_datelife_result()].
#' @inheritParams get_taxon_summary
# get_taxon_summary has param datelife_result
#' @inheritParams patristic_matrix_MRCA
# patristic_matrix_MRCA has param na.rm
#' @return A named numeric vector of MRCA ages for each element given in `datelife_result`.
datelife_result_MRCA <- function(datelife_result, na_rm = TRUE) {
  ages <- sapply(datelife_result, patristic_matrix_MRCA, na_rm = na_rm)
  return(ages)
}

#' Get time of MRCA from patristic matrix. Used in [datelife_result_MRCA()].
#' @param patristic_matrix A patristic matrix (aka a `datelifeResult` object of length 1)
#' @param na_rm If `TRUE`, it drops rows containing `NA`s from the `datelifeResult`
#'   patristic matrix; if `FALSE`, it returns `NA` where there are missing entries.
#' @return The depth of the MRCA as a numeric vector.
patristic_matrix_MRCA <- function(patristic_matrix, na_rm = TRUE) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic_matrix, na.rm = na_rm))
}

datelife_result_newick_all <- function(datelife_result, na_rm = TRUE){
	trees <- sapply(datelife_result, patristic_matrix_to_newick)
	if (na_rm) {
		trees <- trees[which(!is.na(trees))]
	}
	return(trees)
}

datelife_result_phylo_all <- function(datelife_result, na_rm = TRUE){
	trees <- suppressWarnings(lapply(datelife_result, patristic_matrix_to_phylo))
	# suppress warning "Converting from patristic distance matrix to a tree resulted in some negative branch lengths"
	if (na_rm) {
		trees <- trees[which(!is.na(trees))]
	}
	if (length(trees) == 0 ) {
		message("Failed to summarize to 'phylo_all'")
		return(NA)
	}
	if (length(trees) == 1) {
		class(trees) <- "phylo"
	} else {
		class(trees) <- "multiPhylo"
	}
	return(trees)
}

datelife_result_median_and_sdm <- function(datelife_result,
	                                         datelife_query,
																					 criterion = "taxa",
																					 na_rm = TRUE) {
	best_grove <- get_best_grove(datelife_result, criterion = criterion, n = 2)$best_grove
	if (inherits(datelife_query, "datelifeQuery")) {
    target_tree <- datelife_query$phy
  } else {
    target_tree <- NULL
  }
  median_phylo <- datelife_result_median(best_grove, target_tree = target_tree)
	median_newick <- ape::write.tree(median_phylo)
  sdm_phylo <- datelife_result_sdm_phylo(best_grove, target_tree = target_tree)
	sdm_newick <- ape::write.tree(sdm_phylo)
  return(list(median_phylo = median_phylo,
		          median_newick = median_newick,
							sdm_phylo = sdm_phylo,
							sdm_newick = sdm_newick))
}

datelife_result_html <- function(datelife_result, datelife_query, na_rm = TRUE) {
	taxon_summ <- get_taxon_summary(
    datelife_result = datelife_result,
    datelife_query = datelife_query
  )
	out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
	out.vector1 <- paste(out.vector1, "</th></tr>", sep = "")
	ages <- datelife_result_MRCA(datelife_result, na_rm = na_rm)
	trees <- sapply(datelife_result, patristic_matrix_to_newick)
	out.vector2 <- c()
	for (result.index in sequence(length(datelife_result))) {
		out.vector2 <- paste(out.vector2, "<tr><td>", ages[result.index], "</td><td>", sum(!is.na(diag(datelife_result[[result.index]]))), "</td><td>", names(datelife_result)[result.index], "</td><td>", trees[result.index], sep = "")
		out.vector2 <- paste(out.vector2, "</td></tr>", sep = "")
	}
	out.vector <- paste(out.vector1, out.vector2, "</table>", sep = "")
		out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
		for (i in 1:length(taxon_summ$absent_taxa)) {
			out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", taxon_summ$absent_taxa[i], "</td><tr>", sep = "")
		}
		out.vector4 <- paste(out.vector4, "</table>", sep = "")
		out.vector <- paste(out.vector, out.vector4, sep = "")
	return.object <- out.vector
}

datelife_result_data_frame <-  function(datelife_result, na_rm = TRUE){
	out.df <- data.frame()
	ages <- datelife_result_MRCA(datelife_result, na_rm = na_rm)
	trees <- sapply(datelife_result, patristic_matrix_to_newick)
	for (result.index in sequence(length(datelife_result))) {
		out.line <- data.frame(Age = ages[result.index],
			                     Ntax = sum(!is.na(diag(datelife_result[[result.index]]))),
													 Citation = names(datelife_result)[result.index],
													 Newick = trees[result.index])
		if (result.index == 1) {
			out.df <- out.line
		} else {
			out.df <- rbind(out.df, out.line)
		}
	}
	rownames(out.df) <- NULL
	return.object <- out.df
}
