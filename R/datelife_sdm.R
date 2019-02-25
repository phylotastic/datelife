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
			phy <- sdm_matrix_to_phylo(SDM.result)$phy # this also contains the age distributions and calibrations used
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
#' Go from an SDM matrix to anultrametric phylo object.
#' @param sdm_matrix A matrix from get_sdm
#' @return An ultrametric phylo object.
#' @export
sdm_matrix_to_phylo <- function(sdm_matrix){ # enhance: allow other methods, not only bladj.
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
	sdm_matrix <- sdm_matrix*0.5  # bc it's total distance tip to tip
	ages <- tA <- tB <- c() # compute the final length of the data frame: it's ncol(xx)^2 - sum(1:(ncol(xx)-1))
	# calibrations <- matrix(nrow = ncol(xx)^2 - sum(1:(ncol(xx)-1)), ncol = 3)
	# start <- ?
	# identify if SDM matrix has some negative values; extract taxon names:
	negs <- which(sdm_matrix < 0)
	neg_names <- rownames(sdm_matrix)[ceiling(negs/nrow(sdm_matrix))]
	# extract unique ages from sdm_matrix:
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
	# enhance: where does this negative values come from in SDM?
	# get a backbone tree:
	# chronogram <- geiger::PATHd8.phylo(phy_target, calibrations)
	# try(chronogram <- geiger::PATHd8.phylo(phy_target, calibrations), silent = TRUE)
	target_tree <- patristic_matrix_to_phylo(sdm_matrix)
	target_tree$edge.length <- NULL
	target_tree <- tree_add_nodelabels(tree = target_tree, node_index = "consecutive")  # all nodes need to be named so make_bladj_tree runs properly
	# get the coincident nodes:
	target_tree_nodes <- sapply(seq(nrow(calibrations)), function(i)
			phytools::findMRCA(tree = target_tree,
								 tips = as.character(calibrations[i,c("taxonA", "taxonB")]),
								 type = "node"))
	target_tree_nodes <- target_tree_nodes - ape::Ntip(target_tree)
	all_nodes <- sort(unique(target_tree_nodes))
	# get the node age distribution:
	all_ages <- lapply(all_nodes, function(i) calibrations[target_tree_nodes == i, "Age"])
	# any(sapply(all_ages, is.null)) # if FALSE, all nodes have at least one calibration.
	calibrations2 <- data.frame(MRCA = paste0("n", all_nodes), MinAge = sapply(all_ages, min), MaxAge= sapply(all_ages, max), node = all_nodes)
	# calibrations2$MRCA is a factor so have to be made as.character to work with bladj
	new.phy <- make_bladj_tree(tree = target_tree, nodenames = as.character(calibrations2$MRCA), nodeages = sapply(seq(nrow(calibrations2)), function(i) sum(calibrations2[i,c("MinAge", "MaxAge")])/2))
	new.phy$clustering_method <- "sdm"
	return(list(phy = new.phy, sdm_ages_distribution = all_ages, calibrations = calibrations2))
}
