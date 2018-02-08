#  datelife internal functions

#' Summarize patristic matrix array (by default, median). Used inside: summarize_datelife_result.
#' @param patristic.matrix.array 3D array of patristic matrices
#' @param fn The function to use to summarize
#' @return A 2d array with the median (or max, or mean, etc) of the input array
#' @export
SummaryPatristicMatrixArray <- function(patristic.matrix.array, fn = stats::median) {
  return(apply(patristic.matrix.array, MARGIN = c(1,2), fn, na.rm = TRUE))
}

#' Convert patristic matrix to a newick string. Used inside: summarize_datelife_result.
#' @param patristic.matrix A patristic matrix
#' @return A newick string
#' @export
patristic_matrix_to_newick <- function(patristic.matrix) {
  tree <- patristic_matrix_to_phylo(patristic.matrix)
  if(class(tree)=="phylo") {
  	return(ape::write.tree(tree))
  }
  return(NA)
}

#' Find the index of relevant studies in a opentree_chronograms object. Used inside: summarize_datelife_result.
#' @param datelife_result The patristic.matrices that will be used
#' @param cache The cache of studies
#' @return A vector with the indices of studies that have relevant info
#' @export
FindMatchingStudyIndex <- function(datelife_result, cache = get("opentree_chronograms")) {
    return(which(names(cache$trees) %in% names(datelife_result)))
}

#' Get time of MRCA from patristic matrix. Used in: datelife_result_MRCA.
#' @param patristic.matrix A patristic matrix
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return The depth of the MRCA
#' @export
patristic_matrix_MRCA <- function(patristic.matrix, partial = TRUE) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic.matrix, na.rm = partial))
}

#' Get vector of MRCAs from a datelifeResult object. Used in: summarize_datelife_result.
#' @param datelife_result An object from get_datelife_result function.
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return Vector of MRCA ages with names same as in datelife_result
#' @export
datelife_result_MRCA <- function(datelife_result, partial = TRUE) {
	ages <- sapply(datelife_result, patristic_matrix_MRCA, partial = partial)
	return(ages)
}

#' Convert patristic matrix to a phylo object. Used inside: summarize_datelife_result, CongruiyTree.
#' @param patristic.matrix A patristic matrix
#' @return A rooted phylo object
#' @export
patristic_matrix_to_phylo <- function(patristic.matrix) {
  if(anyNA(patristic.matrix)) {
  	patristic.matrix <- patristic.matrix[rowSums(is.na(patristic.matrix)) != ncol(patristic.matrix),colSums(is.na(patristic.matrix)) != nrow(patristic.matrix)]
  }
  if(dim(patristic.matrix)[1] < 2) {
  	return(NA)
  }
	tree <- NA
	if(dim(patristic.matrix)[1] == 2) {
		tree <- ape::rtree(n = 2, rooted = TRUE, tip.label = rownames(patristic.matrix), br = patristic.matrix[1,2]/2)
	} else {
  	tree <- ape::nj(patristic.matrix)
	}
  if(ape::Ntip(tree)>2) {
    tree <- 	phangorn::midpoint(tree)
  }
	if(length(which(tree$edge.length<0))>0) {
		warning(paste("Converting from patristic distance matrix to a tree resulted in some negative branch lengths; the largest by magnitude was", min(tree$edge.length)))
		tree$edge.length[which(tree$edge.length<0)] <- 0 #sometimes NJ returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
	}
  return(tree)
}

#' Figure out which subset function to use. Used inside: get_datelife_result
#' @param study.element The thing being passed in: an array or a phylo to serve as reference
#' @param taxa Vector of taxon names to get a subset for
#' @param phy A user tree to congruify in phylo format (ape)
#' @param phy4 A user tree to congruify in phylo4 format (phylobase)
#' @inheritParams datelife_search
#' @return A patristic matrix with for the taxa.
#' @export
GetSubsetArrayDispatch <- function(study.element, taxa, phy = NULL, phy4 = NULL, dating_method = "PATHd8") {
  if(class(study.element) == "array") {
    return(GetSubsetArrayBoth(study.element, taxa, phy, phy4, dating_method))
  } else {
    return(GetSubsetArrayBothFromPhylo(reference.tree.in = study.element, taxa.in = taxa, phy.in = phy, phy4.in = phy4, dating_method.in = dating_method))
  }
}

#' Take results.list and process it. Used inside: get_datelife_result
#' @param results.list A list returned from using GetSubsetArrayDispatch on opentree_chronograms$trees
#' @param taxa A vector of taxa to match
#' @param partial If TRUE, return matrices that have only partial matches
#' @return A list with the patristic.matrices that are not NA
#' @export
ProcessResultsList <- function(results.list, taxa = NULL, partial = FALSE) {
	if(is.null(taxa)) {
		taxa <- unique(unname(unlist(lapply(final.matrices, rownames))))
	}
	patristic.matrices <- lapply(results.list, "[[", "patristic.matrix.array")

	final.matrices <- patristic.matrices[!is.na(patristic.matrices)]

	if(!partial) {
		final.matrices <- final.matrices[sapply(final.matrices, AllMatching, taxa = taxa)]
	}
	if(length(final.matrices)>0) {
		to.delete <- c()
		for (i in sequence(length(final.matrices))) {
			if(all(is.na(final.matrices[[i]]))) {
				to.delete <- c(to.delete, i)
			}
		}
		if(length(to.delete)>0) {
			final.matrices <- final.matrices[-to.delete]
		}
	}
	return(final.matrices)
}

#' Are all desired taxa in the patristic.matrix? Used inside: ProcessResultsList.
#' @param patristic.matrix A patristic matrix, rownames and colnames must be taxa
#' @param taxa A vector of taxon names
#' @return A Boolean
#' @export
AllMatching <- function(patristic.matrix, taxa) {
	return(sum(!(taxa %in% rownames(patristic.matrix) ))==0)
}

#Used inside: GetSubsetArrayDispatch.
GetSubsetArrayBoth <- function(patristic.matrix.array, taxa, phy = NULL, phy4 = NULL, dating_method = "PATHd8") {
  if (is.null(phy)) {
    return(GetSubsetArray(patristic.matrix.array = patristic.matrix.array, taxa = taxa, phy4 = phy4))
  }
  else { #congruify
    return(GetSubsetArrayCongruify(patristic.matrix.array = patristic.matrix.array, taxa = taxa, phy = phy, dating_method))
  }
}

# Used inside: GetSubsetArrayBoth
GetSubsetArray <- function(patristic.matrix.array, taxa, phy4 = NULL) {
  #gets a subset of the patristic.matrix.array. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa,colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array = patristic.matrix.array,problem = problem))
}

# Used inside: GetSubsetArrayBoth.
GetSubsetArrayCongruify <- function(patristic.matrix.array, taxa, phy = NULL, dating_method = "PATHd8") {
  #gets a subset of the patristic.matrix.array.
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa, colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array = patristic.matrix.array, problem = problem))

    }
  }
  patristic.matrix.list <- SplitArray(patristic.matrix.array)
  patristic.matrix.array <- datelife_result_bind(lapply(patristic.matrix.list, CongruifyTree, query.tree = phy, scale = dating_method)) #yes, this should be parallel
  return(list(patristic.matrix.array = patristic.matrix.array, problem = problem))
}

# Used inside: GetSubsetArrayCongruify
SplitArray <- function(patristic.matrix.array) {
  return(lapply(sequence(dim(patristic.matrix.array)[3]),AsubForLapply,patristic.matrix.array))
}

# Used inside: SplitArray
AsubForLapply <- function(idx, x, dims = 3) {
  return(abind::asub(x, idx, dims))
}

#' Convert list of patristic matrices to a 3D array. Used inside: summarize_datelife_result, GetSubsetArrayCongruify.
#' @param patristic.matrix.list List of patristic matrices
#' @param pad If TRUE, pad missing entries
#' @return A 3d array of patristic matrices
#' @export
datelife_result_bind <- function(patristic.matrix.list, pad = TRUE) {
  all.taxa <- sort(unique(unname(unlist(lapply(patristic.matrix.list, rownames)))))
  if(pad) {
    patristic.matrix.list <- lapply(patristic.matrix.list, datelife_result_pad, all.taxa = all.taxa)
  }
  original.size<-length(patristic.matrix.list)
  patristic.matrix.list<-lapply(patristic.matrix.list,patristic_matrix_reorder)
  if(length(patristic.matrix.list)<1) {
    stop(paste("The patristic matrices you are trying to bind are too few; input was ", original.size, " and current length is ", length(patristic.matrix.list), sep = ""))
  }
  standard.rownames<-rownames(patristic.matrix.list[[1]])
  standard.colnames<-colnames(patristic.matrix.list[[1]])
  matching.names<-sapply(patristic.matrix.list,TestNameOrder,standard.rownames,standard.colnames)
  if (sum(matching.names) != length(matching.names)) {
    stop("The patristic matrices you are trying to bind do not have the same taxa")
  }
  return(abind::abind(patristic.matrix.list, along = 3 ))
}

#' Function to fill in empty cells in a patristic matrix for missing taxa. Used in: datelife_result_bind.
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @param all.taxa A vector of the names of all taxa you want, including ones not in the patristic matrix
#' @return Patristic.matrix for all.taxa, with NA for entries between taxa where at least one was not in the original patristic.matrix
#' @export
datelife_result_pad <- function(patristic.matrix, all.taxa) {
	number.missing <- length(all.taxa) - dim(patristic.matrix)[1]
	final.matrix <- patristic.matrix
	if(number.missing>0) {
		final.matrix <- rbind(patristic.matrix, matrix(nrow = number.missing, ncol = dim(patristic.matrix)[2]))
		final.matrix <- cbind(final.matrix, matrix(ncol = number.missing, nrow = dim(final.matrix)[1]))
		rownames(final.matrix) <- c(rownames(patristic.matrix), all.taxa[-which(all.taxa %in% rownames(patristic.matrix))])
	 	colnames(final.matrix) <- c(colnames(patristic.matrix), all.taxa[-which(all.taxa %in% colnames(patristic.matrix))])
	}
 	return(patristic_matrix_reorder(final.matrix))
}

#' Function to reorder a matrix so that row and column labels are in alphabetical order. Used in: datelife_result_pad.
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return patristic.matrix A patristic matrix with row and column names for taxa in alphabetial order
#' @export
patristic_matrix_reorder <- function(patristic.matrix) {
  return(patristic.matrix[order(rownames(patristic.matrix)),order(colnames(patristic.matrix))])
}

# Used in: datelife_result_bind.
TestNameOrder <- function(patristic.matrix, standard.rownames, standard.colnames) {
  if (compare::compare(rownames(patristic.matrix),standard.rownames)$result!= TRUE) {
    return(FALSE)
  }
  if (compare::compare(colnames(patristic.matrix),standard.colnames)$result!= TRUE) {
    return(FALSE)
  }
  return(TRUE)
}

# Used inside: GetSubsetArrayCongruify.
CongruifyTree <- function(patristic.matrix, query.tree, dating_method = "PATHd8", attempt.fix = TRUE) {
  	result.matrix <- matrix(nrow = dim(patristic.matrix)[1], ncol = dim(patristic.matrix)[2])
  	if(is.null(query.tree$edge.length)) {
    	query.tree$edge.length<-numeric(nrow(query.tree$edge))
  	}
#	try(result.matrix<-phylo_to_patristic_matrix(phylo_tiplabel_underscore_to_space(geiger::congruify.phylo(phylo_tiplabel_space_to_underscore(patristic_matrix_to_phylo(patristic.matrix)), phylo_tiplabel_space_to_underscore(query.tree), NULL, 0, scale = dating_method)$phy)))
	try(result.matrix <- phylo_to_patristic_matrix(CongruifyAndCheck(reference = patristic_matrix_to_phylo(patristic.matrix), target = query.tree, scale = dating_method, attempt.fix = attempt.fix)))
  return(result.matrix)
}

#Note that originally trees were stored as patristic matrices. This was intended
#to make subsetting fast. The downside is large memory usage. Klaus Schliep wrote
#fast tree subsetting for phylo and multiphylo objects, so now trees are stored
#internally as objects of this type, but with the final output after pruning
#going through patristic matrices.

#Some trees are so large that they can't be stored as patristic distance matrices. For all others,
#patristic matrices are better. For example, for the 20,000 HeathEtAl2012 trees of 35 taxa,
#getting a subset down to two taxa takes 0.0475 seconds just for the pruning, 0.0504 seconds
#for pruning and getting a subset, for a single tree (run times go up linearly with number of trees:
#pruning and converting 1000 trees takes 3 seconds). Subsetting 1000 trees from the patristic
#distance matrix takes just 0.0013 seconds.

#in case we want to cache. Not clear we do.
# Used inside: CongruifyTree, GetSubsetArrayFromPhylo and CongruifyTreeFromPhylo
phylo_to_patristic_matrix <- function(phy, test = TRUE, tol = 0.01, option = 2) {
	# stores the distance between taxa
	patristic.matrix <- NA
	if(class(phy) == "phylo") {
		if (test) {
			if (!ape::is.ultrametric(phy, tol = tol, option = option)) {
				stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
			}
		}
		patristic.matrix <- stats::cophenetic(phy)
	}
	return(patristic.matrix)
}

#Used inside: GetSubsetArrayDispatch.
GetSubsetArrayBothFromPhylo <- function(reference.tree.in, taxa.in, phy.in = NULL, phy4.in = NULL, dating_method.in = "PATHd8") {
#COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
#  if (class(reference.tree)=="phylo") {
#    reference.tree<-c(reference.tree) #from here in, assumes multiphylo object, even if a single tree
#  }
	congruify = FALSE
	if(!is.null(phy.in[1])) {
		congruify = TRUE
		if(is.na(phy.in[1])) {
			congruify = FALSE
		}
	}
  if (!congruify) {
    return(GetSubsetArrayFromPhylo(reference.tree = reference.tree.in, taxa = taxa.in, phy4 = phy4.in, dating_method.in))
  }
  else { #congruify
    return(GetSubsetArrayCongruifyFromPhylo(reference.tree = reference.tree.in, taxa = taxa.in, phy = phy.in, dating_method.in))
  }

}

# Used inside: GetSubsetArrayBothFromPhylo.
GetSubsetArrayFromPhylo <- function(reference.tree, taxa, phy4 = NULL, dating_method="PATHd8") {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
    reference.tree <- PruneTree(reference.tree, taxa)
    #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem <- "none"
  patristic.matrix.array <- NA
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
  	patristic.matrix.array <- phylo_to_patristic_matrix(reference.tree)
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array = patristic.matrix.array,problem = problem))
}

# Used inside: GetSubsetArrayBothFromPhylo.
GetSubsetArrayCongruifyFromPhylo <- function(reference.tree, taxa, phy = NULL, dating_method = "PATHd8") {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
   reference.tree<-PruneTree(reference.tree, taxa)
   #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem.new <- "none"
  patristic.matrix.array.new <- NA
  if (final.size < length(taxa)) {
    problem.new <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array.new <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array = patristic.matrix.array.new, problem = problem.new))

    }
  }
  if (final.size >= 3) {
  	patristic.matrix.array.new <- CongruifyTreeFromPhylo(reference.tree, query.tree = phy, dating_method = dating_method)
  }
  return(list(patristic.matrix.array= patristic.matrix.array.new,problem= problem.new))
}

# Used inside: GetSubsetArrayFromPhylo and GetSubsetArrayCongruifyFromPhylo.
PruneTree <- function(phy, taxa) {
	return(ape::drop.tip(phy, tip = phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

# Used inside: GetSubsetArrayCongruifyFromPhylo.
CongruifyTreeFromPhylo <- function(reference.tree, query.tree, dating_method = "PATHd8", attempt.fix = TRUE) {
  result.matrix <- matrix(nrow = ape::Ntip(reference.tree), ncol = ape::Ntip(reference.tree))
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
	try(result.matrix <- phylo_to_patristic_matrix(CongruifyAndCheck(reference = reference.tree, target = query.tree, scale = dating_method, attempt.fix = attempt.fix)))
  return(result.matrix)
}

# Used inside: CongruifyTree and CongruifyTreeFromPhylo.
CongruifyAndCheck <- function(reference, target, taxonomy = NULL, tol = 0.01, option = 2, scale = "pathd8", attempt.fix = TRUE) {
  if(!ape::is.ultrametric(reference, tol = tol, option = option)) {
    return(NA)
  }
	new.tree <- phylo_tiplabel_underscore_to_space(suppressWarnings(geiger::congruify.phylo(phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(target), taxonomy = taxonomy, tol = tol, scale = scale)$phy)) #suppressing warnings b/c geiger ignores tolerance
	if(anyNA(new.tree$edge.length) & attempt.fix) {
		warning("Congruification resulted in NA edge lengths. Resolving polytomies and making up starting branch lengths")
		new.tree <- phylo_tiplabel_underscore_to_space(geiger::congruify.phylo(phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(ape::compute.brlen(ape::multi2di(target))), taxonomy, tol, scale)$phy)
		if(anyNA(new.tree$edge.length)) {
			new.tree <- NA
		}
	}
	new.tree$edge.length[which(new.tree$edge.length<0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
	return(new.tree)
}


#' Convert spaces to underscores in trees. Used in: make_mrbayes_runfile, get_mrbayes_node_calibrations, phylo_get_singleton_outgroup, CongruifyAndCheck, CongruifyTree.
#' @param phy A phylo object
#' @return A phylo object
#' @export
phylo_tiplabel_space_to_underscore <- function(phy) {
	phy$tip.label <- gsub(" ", "_", phy$tip.label)
	return(phy)
}

#' Convert underscores to spaces in trees. Used inside: CongruifyTree, CongruifyAndCheck.
#' @param phy A phylo object
#' @return A phylo object
#' @export
phylo_tiplabel_underscore_to_space <- function(phy) {
	phy$tip.label <- gsub("_", " ", phy$tip.label)
	return(phy)
}
#' Function to remove missing taxa from a datelifeResult object. Used in: datelife_result_sdm.
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return Patristic.matrix for all.taxa
#' @export
patristic_matrix_unpad <- function(patristic.matrix) {
	bad.ones <- which(apply(is.na(patristic.matrix),2,all))
	if(length(bad.ones)>0) {
		patristic.matrix <- patristic.matrix[-bad.ones, -bad.ones]
	}
	return(patristic.matrix)
}
