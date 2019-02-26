#' Function to compute the SDM supertree (Criscuolo et al. 2006) from a datelifeResult object.
#' @inheritParams datelife_result_check
#' @inheritParams get_sdm
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
	# add check datelife_result
	phy <- NA
	# used.studies <- names(datelife_result)
	if(length(datelife_result) == 1){
		phy <- patristic_matrix_to_phylo(datelife_result, fix_negative_brlen = FALSE)
		unpadded.matrices <- datelife_result
	} else {
		# datelife_result <- best_grove
		unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
		good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose)
		if(length(good.matrix.indices) > 1) {
			unpadded.matrices <- unpadded.matrices[good.matrix.indices]
			SDM.result <- get_sdm(unpadded.matrices, weighting, verbose)
			# it is important to use upgma as clustering method; nj produces much younger ages when the matrix comes from sdm
			# phy <- patristic_matrix_to_phylo(SDM.result, clustering_method = clustering_method, fix_negative_brlen = TRUE)
			phy <- summary_matrix_to_phylo(SDM.result)$phy # this also contains the age distributions and calibrations used
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
#' Get SDM result for list of working matrices.
#' @param unpadded.matrices A list of patristic matrices, a datelifeResult object.
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
#' @inheritParams datelife_search
#' @return A matrix.
#' @export
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
#' Get good matrices for SDM
#' @inheritParams get_sdm
#' @return A numeric vector of good matrix indices in unpadded.matrices.
#' @export
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
