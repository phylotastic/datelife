#  datelife internal functions

#' Summarize patristic matrix array (by default, median). Used inside: summarize_datelife_result.
#' @param patristic_matrix_array 3D array of patristic matrices
#' @param fn The function to use to summarize
#' @return A 2d array with the median (or max, or mean, etc) of the input array
summary_patristic_matrix_array <- function(patristic_matrix_array, fn = stats::median) {
  return(apply(patristic_matrix_array, MARGIN = c(1,2), fn, na.rm = TRUE))
}



#' Find the index of relevant studies in a opentree_chronograms object. Used inside: summarize_datelife_result.
#' @inheritParams datelife_result_check
#' @param cache The cache of studies
#' @return A vector with the indices of studies that have relevant info
datelife_result_study_index <- function(datelife_result,
                                        cache = "opentree_chronograms") {

  if("opentree_chronograms" %in% cache){
   utils::data("opentree_chronograms")
   cache <- get("opentree_chronograms")
  }
  return(which(names(cache$trees) %in% names(datelife_result)))
}

#' Get time of MRCA from patristic matrix. Used in: datelife_result_MRCA.
#' @param patristic_matrix A patristic matrix
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return The depth of the MRCA
patristic_matrix_MRCA <- function(patristic_matrix, partial = TRUE) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic_matrix, na.rm = partial))
}

#' Get vector of MRCAs from a datelifeResult object. Used in: summarize_datelife_result.
#' @inheritParams datelife_result_check
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return Vector of MRCA ages with names same as in datelife_result
datelife_result_MRCA <- function(datelife_result, partial = TRUE) {
	ages <- sapply(datelife_result, patristic_matrix_MRCA, partial = partial)
	return(ages)
}

#' Figure out which subset function to use. Used inside: get_datelife_result
#' @param study_element The thing being passed in: an array or a phylo object to serve as reference for congruification
#' @param taxa Vector of taxon names to get a subset for
#' @param phy A user tree to congruify as phylo object (ape)
#' @param phy4 A user tree to congruify in phylo4 format (phylobase)
#' @param dating_method The method used for tree dating.
#' @return A patristic matrix with ages for the target taxa.
get_subset_array_dispatch <- function(study_element, taxa, phy = NULL, phy4 = NULL, dating_method = "PATHd8") {
  if(class(study_element) == "array") {
    return(patristic_matrix_array_subset_both(study_element, taxa, phy, phy4, dating_method))
  } else {
    return(phylo_subset_both(reference_tree.in = study_element, taxa.in = taxa, phy.in = phy, phy4.in = phy4, dating_method.in = dating_method))
  }
}

#' Take results_list and process it. Used inside: get_datelife_result
#' @param results_list A list returned from using get_subset_array_dispatch on opentree_chronograms$trees
#' @param taxa A vector of taxa to match
#' @param partial If TRUE, return matrices that have only partial matches
#' @return A list with the patristic.matrices that are not NA
results_list_process <- function(results_list, taxa = NULL, partial = FALSE) {
	if(is.null(taxa)) {
		taxa <- unique(unname(unlist(lapply(final_matrices, rownames))))
	}
	patristic.matrices <- lapply(results_list, "[[", "patristic_matrix_array")

	final_matrices <- patristic.matrices[!is.na(patristic.matrices)]

	if(length(final_matrices) > 0) {
        if(!partial) {
    		final_matrices <- final_matrices[sapply(final_matrices, patristic_matrix_taxa_all_matching, taxa = taxa)]
    	}
		to.delete <- c()
		for (i in sequence(length(final_matrices))) {
			if(all(is.na(final_matrices[[i]]))) {
				to.delete <- c(to.delete, i)
			}
		}
		if(length(to.delete) > 0) {
			final_matrices <- final_matrices[-to.delete]
		}
	}
	return(final_matrices)
}

#' Are all desired taxa in the patristic_matrix? Used inside: results_list_process.
#' @param patristic_matrix A patristic matrix, rownames and colnames must be taxa
#' @param taxa A vector of taxon names
#' @return A Boolean
patristic_matrix_taxa_all_matching <- function(patristic_matrix, taxa) {
	return(sum(!(taxa %in% rownames(patristic_matrix) )) == 0)
}

#Used inside: get_subset_array_dispatch.
patristic_matrix_array_subset_both <- function(patristic_matrix_array, taxa, phy = NULL, phy4 = NULL, dating_method = "PATHd8") {
  if (is.null(phy)) {
    return(patristic_matrix_array_subset(patristic_matrix_array = patristic_matrix_array, taxa = taxa, phy4 = phy4))
  }
  else { #congruify
    return(patristic_matrix_array_congruify(patristic_matrix_array = patristic_matrix_array, taxa = taxa, phy = phy, dating_method))
  }
}

# Used inside: patristic_matrix_array_subset_both
patristic_matrix_array_subset <- function(patristic_matrix_array, taxa, phy4 = NULL) {
  #gets a subset of the patristic_matrix_array. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic_matrix_array <- patristic_matrix_array[rownames(patristic_matrix_array) %in% taxa, colnames(patristic_matrix_array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic_matrix_array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic_matrix_array = patristic_matrix_array,problem = problem))
}

# Used inside: patristic_matrix_array_subset_both. patristic_matrix_array_congruify
patristic_matrix_array_congruify <- function(patristic_matrix_array, taxa, phy = NULL,
    dating_method = "PATHd8") {
  #gets a subset of the patristic_matrix_array.
  patristic_matrix_array <- patristic_matrix_array[rownames(patristic_matrix_array) %in% taxa, colnames(patristic_matrix_array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic_matrix_array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))

    }
  }
  patristic_matrix_list <- patristic_matrix_array_split(patristic_matrix_array)
  patristic_matrix_array <- patristic_matrix_list_to_array(lapply(patristic_matrix_list, patristic_matrix_array_phylo_congruify, target_tree = phy, scale = dating_method)) #yes, this should be parallel
  return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))
}

# Used inside: patristic_matrix_array_congruify
patristic_matrix_array_split <- function(patristic_matrix_array) {
  return(lapply(sequence(dim(patristic_matrix_array)[3]),asub_for_lapply,patristic_matrix_array))
}

# Used inside: patristic_matrix_array_split to asub.list?
asub_for_lapply <- function(idx, x, dims = 3) {
  return(abind::asub(x, idx, dims))
}

#' Convert list of patristic matrices to a 3D array. Used inside: summarize_datelife_result, patristic_matrix_array_congruify.
#' @param patristic_matrix_list List of patristic matrices
#' @param pad If TRUE, pad missing entries
#' @return A 3d array of patristic matrices
patristic_matrix_list_to_array <- function(patristic_matrix_list, pad = TRUE) {
  all_taxa <- sort(unique(unname(unlist(lapply(patristic_matrix_list, rownames)))))
  if(pad) {
    patristic_matrix_list <- lapply(patristic_matrix_list, patristic_matrix_pad, all_taxa = all_taxa)
  }
  original.size <- length(patristic_matrix_list)
  patristic_matrix_list <- lapply(patristic_matrix_list,patristic_matrix_name_reorder)
  if(length(patristic_matrix_list) < 1) {
    stop(paste0("The patristic matrices you are trying to bind are too few; input was ", original.size, " and current length is ", length(patristic_matrix_list)))
  }
  standard.rownames <- rownames(patristic_matrix_list[[1]])
  standard.colnames <- colnames(patristic_matrix_list[[1]])
  matching.names <- sapply(patristic_matrix_list,patristic_matrix_name_order_test,standard.rownames,standard.colnames)
  if (sum(matching.names) != length(matching.names)) {
    stop("The patristic matrices you are trying to bind do not have the same taxa")
  }
  return(abind::abind(patristic_matrix_list, along = 3 ))
}

#' Function to fill in empty cells in a patristic matrix for missing taxa. Used in: patristic_matrix_list_to_array.
#' @param patristic_matrix A patristic matrix with row and column names for taxa
#' @param all_taxa A vector of the names of all taxa you want, including ones not in the patristic matrix
#' @return patristic_matrix for all_taxa, with NA for entries between taxa where at least one was not in the original patristic_matrix
patristic_matrix_pad <- function(patristic_matrix, all_taxa) {
	number.missing <- length(all_taxa) - dim(patristic_matrix)[1]
	final_matrix <- patristic_matrix
	if(number.missing>0) {
		final_matrix <- rbind(patristic_matrix, matrix(nrow = number.missing, ncol = dim(patristic_matrix)[2]))
		final_matrix <- cbind(final_matrix, matrix(ncol = number.missing, nrow = dim(final_matrix)[1]))
		rownames(final_matrix) <- c(rownames(patristic_matrix), all_taxa[-which(all_taxa %in% rownames(patristic_matrix))])
	 	colnames(final_matrix) <- c(colnames(patristic_matrix), all_taxa[-which(all_taxa %in% colnames(patristic_matrix))])
	}
 	return(patristic_matrix_name_reorder(final_matrix))
}

#' Function to reorder a matrix so that row and column labels are in alphabetical order. Used in: patristic_matrix_pad.
#' @param patristic_matrix A patristic matrix with row and column names for taxa
#' @return patristic_matrix A patristic matrix with row and column names for taxa in alphabetical order
patristic_matrix_name_reorder <- function(patristic_matrix) {
  return(patristic_matrix[order(rownames(patristic_matrix)),order(colnames(patristic_matrix))])
}

# Used in: patristic_matrix_list_to_array.
patristic_matrix_name_order_test <- function(patristic_matrix, standard.rownames,
    standard.colnames) {

  if (compare::compare(rownames(patristic_matrix),standard.rownames)$result!= TRUE) {
    return(FALSE)
  }
  if (compare::compare(colnames(patristic_matrix),standard.colnames)$result!= TRUE) {
    return(FALSE)
  }
  return(TRUE)
}

# Used inside: patristic_matrix_array_congruify.
patristic_matrix_array_phylo_congruify <- function(patristic_matrix, target_tree,
    dating_method = "PATHd8", attempt.fix = TRUE) {

  	result_matrix <- matrix(nrow = dim(patristic_matrix)[1], ncol = dim(patristic_matrix)[2])
  	if(is.null(target_tree$edge.length)) {
    	target_tree$edge.length<-numeric(nrow(target_tree$edge))
  	}
	try(result_matrix <- phylo_to_patristic_matrix(congruify_and_check(
        reference = patristic_matrix_to_phylo(patristic_matrix),
        target = target_tree, scale = dating_method, attempt.fix = attempt.fix)))
  return(result_matrix)
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
# Used inside: patristic_matrix_array_phylo_congruify, phylo_get_subset_array and phylo_congruify
phylo_to_patristic_matrix <- function(phy, test = TRUE, tol = 0.01, option = 2) {
	# stores the distance between taxa
	patristic_matrix <- NA
	if(class(phy) == "phylo") {
		if (test) {
			if (!ape::is.ultrametric(phy, tol = tol, option = option)) {
				stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
			}
		}
		patristic_matrix <- stats::cophenetic(phy)
	}
	return(patristic_matrix)
}

#Used inside: get_subset_array_dispatch.
phylo_subset_both <- function(reference_tree.in, taxa.in, phy.in = NULL, phy4.in = NULL, dating_method.in = "PATHd8") {
#COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
#  if (inherits(reference_tree, "phylo")) {
#    reference_tree<-c(reference_tree) #from here in, assumes multiphylo object, even if a single tree
#  }
	congruify = FALSE
	if(inherits(phy.in, "phylo")) {
            congruify = TRUE
	} else {
        congruify = FALSE
    }
  if (congruify) {
      return(phylo_get_subset_array_congruify(reference_tree = reference_tree.in, taxa = taxa.in, phy = phy.in, dating_method = dating_method.in))
  } else {  # when congruify is FALSE:
      return(phylo_get_subset_array(reference_tree = reference_tree.in, taxa = taxa.in, phy4 = phy4.in, dating_method = dating_method.in))
  }
}

# Used inside: phylo_subset_both, when we don't congruify
phylo_get_subset_array <- function(reference_tree, taxa, phy4 = NULL, dating_method = "PATHd8") {
  final.size <- sum(reference_tree$tip.label %in% taxa) # returns number of matches
  if(final.size >= 2) { #it's worth doing the pruning
    reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa)
    #reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa) #phylo_prune_missing_taxa (PruneTree before) is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem <- "none"
  patristic_matrix_array <- NA
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
  	patristic_matrix_array <- phylo_to_patristic_matrix(reference_tree)
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic_matrix_array = patristic_matrix_array,problem = problem))
}

# Used inside: phylo_subset_both.
phylo_get_subset_array_congruify <- function(reference_tree, taxa, phy = NULL, dating_method = "PATHd8") {
  final.size <- sum(reference_tree$tip.label %in% taxa) # returns number of matches
  if(final.size >= 2) {  # it's worth doing the pruning
   reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa)
   #reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa)  # phylo_prune_missing_taxa (used to be names PruneTree) is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem.new <- "none"
  patristic_matrix_array.new <- NA
  if (final.size < length(taxa)) {
    problem.new <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array.new <- NA # to make sure no one uses the zero by mistake
      return(list(patristic_matrix_array = patristic_matrix_array.new, problem = problem.new))
    }
  }
  if (final.size >= 3) {
  	patristic_matrix_array.new <- phylo_congruify(reference_tree, target_tree = phy, dating_method = dating_method)
  }
  return(list(patristic_matrix_array = patristic_matrix_array.new, problem = problem.new))
}

# Used inside: phylo_get_subset_array and phylo_get_subset_array_congruify.
phylo_prune_missing_taxa <- function(phy, taxa) {
	return(ape::drop.tip(phy, tip = phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

# Used inside: phylo_get_subset_array_congruify.
phylo_congruify <- function(reference_tree, target_tree, dating_method = "PATHd8", attempt.fix = TRUE) {
  result_matrix <- matrix(nrow = ape::Ntip(reference_tree), ncol = ape::Ntip(reference_tree))
  if(is.null(target_tree$edge.length)) {
    target_tree$edge.length <- numeric(nrow(target_tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
	try(result_matrix <- phylo_to_patristic_matrix(congruify_and_check(reference = reference_tree, target = target_tree, scale = dating_method, attempt.fix = attempt.fix)))
  return(result_matrix)
}

# Used inside: patristic_matrix_array_phylo_congruify and phylo_congruify.
congruify_and_check <- function(reference, target, taxonomy = NULL, tol = 0.01,
    option = 2, scale = "pathd8", attempt.fix = TRUE) {
  if(!ape::is.ultrametric(reference, tol = tol, option = option)) {
    return(NA)
  }
  new.tree <- phylo_tiplabel_underscore_to_space(suppressWarnings(geiger::congruify.phylo(
    phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(target),
    taxonomy = taxonomy, tol = tol, scale = scale)$phy)) #suppressing warnings b/c geiger ignores tolerance
  if(anyNA(new.tree$edge.length) & attempt.fix) {
	message("Congruification resulted in NA edge lengths. Resolving polytomies and making up starting branch lengths")
	new.tree <- phylo_tiplabel_underscore_to_space(geiger::congruify.phylo(
        phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(
            ape::compute.brlen(ape::multi2di(target))), taxonomy, tol, scale, ncores = 1)$phy)
	if(anyNA(new.tree$edge.length)) {
        message("There are still NAs in edge lengths; returning NA")
		new.tree <- NA
	}
  }
  new.tree$edge.length[which(new.tree$edge.length < 0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
  return(new.tree)
}


#' Convert spaces to underscores in trees. Used in: make_mrbayes_runfile, get_mrbayes_node_calibrations, tree_get_singleton_outgroup, congruify_and_check, patristic_matrix_array_phylo_congruify.
#' @inheritParams phylo_check
#' @return A phylo object
phylo_tiplabel_space_to_underscore <- function(phy) {
	phy$tip.label <- gsub(" ", "_", phy$tip.label)
	return(phy)
}

#' Convert underscores to spaces in trees. Used inside: patristic_matrix_array_phylo_congruify, congruify_and_check.
#' @inheritParams phylo_check
#' @return A phylo object
phylo_tiplabel_underscore_to_space <- function(phy) {
    # a better name for this function would be underscore_to_blank
    # add method .phylo
    # change tip and node labels
	phy$tip.label <- gsub("_", " ", phy$tip.label)
    # make sure there is only one consecutive blank at a time
	return(phy)
}
