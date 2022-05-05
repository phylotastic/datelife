#  datelife internal functions

#' Summarize patristic matrix array (by default, median). Used inside: summarize_datelife_result.
#' @param patristic_matrix_array 3D array of patristic matrices
#' @param fn The function to use to summarize
#' @return A 2d array with the median (or max, or mean, etc) of the input array
summary_patristic_matrix_array <- function(patristic_matrix_array, fn = stats::median) {
  return(apply(patristic_matrix_array, MARGIN = c(1, 2), fn, na.rm = TRUE))
}



#' Find the index of relevant studies in a cached chronogram database.
#'
#' `datelife_result_study_index` is used in [summarize_datelife_result()].
#'
#' @inheritParams summarize_datelife_result
#' @inheritParams datelife_authors_tabulate
#' @return A vector of indices of studies that have relevant information.
datelife_result_study_index <- function(datelife_result,
                                        cache = "opentree_chronograms") {
  if ("opentree_chronograms" %in% cache) {
    utils::data("opentree_chronograms")
    cache <- get("opentree_chronograms")
  }
  return(which(names(cache$trees) %in% names(datelife_result)))
}

#' Figure out which subset function to use.
#'
#' `get_subset_array_dispatch` is used inside [get_datelife_result()]
#'
#' @param study_element The thing being passed in: an `array` or a `phylo` object
#'  to serve as reference for congruification.
#' @param taxa Vector of taxon names to get a subset for.
#' @param phy A user tree to congruify as `phylo` object (ape).
#' @param phy4 A user tree to congruify in `phylo4` format (phylobase).
#' @param dating_method The method used for tree dating.
#' @return A patristic matrix with ages for the target taxa.
get_subset_array_dispatch <- function(study_element,
                                      taxa, phy = NULL,
                                      phy4 = NULL,
                                      dating_method = "PATHd8") {
  if (inherits(study_element, "array")) {
    return(patristic_matrix_array_subset_both(study_element, taxa, phy, phy4, dating_method))
  } else {
    return(phylo_subset_both(reference_tree = study_element, taxa = taxa, phy = phy, phy4 = phy4, dating_method = dating_method))
  }
}

#' Take results_list and process it.
#'
#' `results_list_process` is used inside: [get_datelife_result()]
#'
#' @param results_list A `list` returned from using [get_subset_array_dispatch()] on `opentree_chronograms$trees`
#' @inheritParams get_subset_array_dispatch
#' @param partial If `TRUE`, return matrices that have only partial matches.
#' @return A list with the patristic.matrices that are not `NA`.
results_list_process <- function(results_list, taxa = NULL, partial = FALSE) {
  if (is.null(taxa)) {
    taxa <- unique(unname(unlist(lapply(final_matrices, rownames))))
  }
  patristic.matrices <- lapply(results_list, "[[", "patristic_matrix_array")

  final_matrices <- patristic.matrices[!is.na(patristic.matrices)]

  if (length(final_matrices) > 0) {
    if (!partial) {
      final_matrices <- final_matrices[sapply(final_matrices, patristic_matrix_taxa_all_matching, taxa = taxa)]
    }
    to.delete <- c()
    for (i in sequence(length(final_matrices))) {
      if (all(is.na(final_matrices[[i]]))) {
        to.delete <- c(to.delete, i)
      }
    }
    if (length(to.delete) > 0) {
      final_matrices <- final_matrices[-to.delete]
    }
  }
  return(final_matrices)
}

#' Are all desired taxa in the patristic matrix?
#'
#' `patristic_matrix_taxa_all_matching` is used inside: [results_list_process()].
#'
#' @param patristic_matrix A patristic matrix, `rownames` and `colnames` must be taxa.
#' @inheritParams get_subset_array_dispatch
#' @return A Boolean.
patristic_matrix_taxa_all_matching <- function(patristic_matrix, taxa) {
  return(sum(!(taxa %in% rownames(patristic_matrix))) == 0)
}


#' Are all desired taxa in the patristic matrix array?
#'
#' `patristic_matrix_array_subset_both` is used inside [get_subset_array_dispatch()].
#'
#' @param patristic_matrix_array A patristic matrix array, `rownames` and `colnames` must be taxa.
#' @inheritParams get_subset_array_dispatch
#' @inherit get_subset_array_dispatch return
patristic_matrix_array_subset_both <- function(patristic_matrix_array, taxa, phy = NULL, phy4 = NULL, dating_method = "PATHd8") {
  if (is.null(phy)) {
    return(patristic_matrix_array_subset(patristic_matrix_array = patristic_matrix_array, taxa = taxa, phy4 = phy4))
  } else { # congruify
    return(patristic_matrix_array_congruify(patristic_matrix_array = patristic_matrix_array, taxa = taxa, phy = phy, dating_method))
  }
}

#' Subset a patristic matrix array
#' @inheritParams patristic_matrix_array_subset_both
#' @return A list with a patristic matrix array and a `$problem` if any.
patristic_matrix_array_subset <- function(patristic_matrix_array, taxa, phy4 = NULL) {
  # gets a subset of the patristic_matrix_array. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic_matrix_array <- patristic_matrix_array[rownames(patristic_matrix_array) %in% taxa, colnames(patristic_matrix_array) %in% taxa, ]
  problem <- "none"
  final.size <- sum(rownames(patristic_matrix_array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))
}


#' `patristic_matrix_array_congruify` is used for patristic_matrix_array_subset_both and patristic_matrix_array_congruify.
#'
#' @inheritParams patristic_matrix_array_subset_both
#' @inherit patristic_matrix_array_subset_both return
patristic_matrix_array_congruify <- function(patristic_matrix_array,
                                             taxa,
                                             phy = NULL,
                                             dating_method = "PATHd8") {
  # gets a subset of the patristic_matrix_array.
  patristic_matrix_array <- patristic_matrix_array[rownames(patristic_matrix_array) %in% taxa, colnames(patristic_matrix_array) %in% taxa, ]
  problem <- "none"
  final.size <- sum(rownames(patristic_matrix_array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))
    }
  }
  patristic_matrix_list <- patristic_matrix_array_split(patristic_matrix_array)
  patristic_matrix_array <- patristic_matrix_list_to_array(lapply(patristic_matrix_list, patristic_matrix_array_phylo_congruify, target_tree = phy, scale = dating_method)) # yes, this should be parallel
  return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))
}

#' Split a patristic matrix array
#' Used inside: patristic_matrix_array_congruify
#' @inheritParams patristic_matrix_array_congruify
#' @return  A patristic matrix 3d array.
patristic_matrix_array_split <- function(patristic_matrix_array) {
  asub_for_lapply <- function(idx, x, dims = 3) {
    return(abind::asub(x, idx, dims))
  }
  return(lapply(sequence(dim(patristic_matrix_array)[3]), asub_for_lapply, patristic_matrix_array))
}

#' Convert list of patristic matrices to a 3D array.
#'
#' `patristic_matrix_list_to_array` us ised inside [summarize_datelife_result()], [patristic_matrix_array_congruify()].
#'
#' @param patristic_matrix_list List of patristic matrices
#' @param pad If TRUE, pad missing entries
#' @return A 3d array of patristic matrices
patristic_matrix_list_to_array <- function(patristic_matrix_list, pad = TRUE) {
  all_taxa <- sort(unique(unname(unlist(lapply(patristic_matrix_list, rownames)))))
  if (pad) {
    patristic_matrix_list <- lapply(patristic_matrix_list, patristic_matrix_pad, all_taxa = all_taxa)
  }
  original.size <- length(patristic_matrix_list)
  patristic_matrix_list <- lapply(patristic_matrix_list, patristic_matrix_name_reorder)
  if (length(patristic_matrix_list) < 1) {
    stop(paste0("The patristic matrices you are trying to bind are too few; input was ", original.size, " and current length is ", length(patristic_matrix_list)))
  }
  standard.rownames <- rownames(patristic_matrix_list[[1]])
  standard.colnames <- colnames(patristic_matrix_list[[1]])
  matching.names <- sapply(patristic_matrix_list, patristic_matrix_name_order_test, standard.rownames, standard.colnames)
  if (sum(matching.names) != length(matching.names)) {
    stop("The patristic matrices you are trying to bind do not have the same taxa")
  }
  return(abind::abind(patristic_matrix_list, along = 3))
}

#' Fill in empty cells in a patristic matrix for missing taxa.
#'
#' Used in: [patristic_matrix_list_to_array()].
#'
#' @inheritParams patristic_matrix_taxa_all_matching
#' @param all_taxa A vector of names of all taxa you want, including ones not
#' in the patristic matrix.
#' @return A patristic matrix, with `NA` for entries between taxa
#' where at least one was not in the original patristic matrix.
patristic_matrix_pad <- function(patristic_matrix, all_taxa) {
  number.missing <- length(all_taxa) - dim(patristic_matrix)[1]
  final_matrix <- patristic_matrix
  if (number.missing > 0) {
    final_matrix <- rbind(patristic_matrix, matrix(nrow = number.missing, ncol = dim(patristic_matrix)[2]))
    final_matrix <- cbind(final_matrix, matrix(ncol = number.missing, nrow = dim(final_matrix)[1]))
    rownames(final_matrix) <- c(rownames(patristic_matrix), all_taxa[-which(all_taxa %in% rownames(patristic_matrix))])
    colnames(final_matrix) <- c(colnames(patristic_matrix), all_taxa[-which(all_taxa %in% colnames(patristic_matrix))])
  }
  return(patristic_matrix_name_reorder(final_matrix))
}

#' Reorder a matrix so that row and column labels are in alphabetical order.
#'
#' `patristic_matrix_name_reorder` is only used in: [patristic_matrix_pad()].
#'
#' @inheritParams patristic_matrix_taxa_all_matching
#' @return A patristic matrix with row and column names for taxa in alphabetical order.
patristic_matrix_name_reorder <- function(patristic_matrix) {
  return(patristic_matrix[order(rownames(patristic_matrix)), order(colnames(patristic_matrix))])
}


#' Test the name order of a patristic matrix so that row and column labels are in alphabetical order.
#'
#' `patristic_matrix_name_order_test` is only used in [patristic_matrix_list_to_array()].
#'
#' @inheritParams patristic_matrix_taxa_all_matching
#' @param standard.rownames A character vector of row names.
#' @param standard.colnames A character vector of column names.
#' @return Boolean.
patristic_matrix_name_order_test <- function(patristic_matrix,
                                             standard.rownames,
                                             standard.colnames) {
  if (compare::compare(rownames(patristic_matrix), standard.rownames)$result != TRUE) {
    return(FALSE)
  }
  if (compare::compare(colnames(patristic_matrix), standard.colnames)$result != TRUE) {
    return(FALSE)
  }
  return(TRUE)
}

# Used inside: patristic_matrix_array_congruify.
#' Congruify a patristic matrix array from a given `phylo` object.
#' @inheritParams patristic_matrix_taxa_all_matching
# # ' @inheritParams patristic_matrix_array_congruify
# # ' @inheritParams phylo_get_subset_array
#' @inheritParams phylo_congruify
#' @inherit phylo_congruify return
patristic_matrix_array_phylo_congruify <- function(patristic_matrix,
                                                   target_tree,
                                                   dating_method = "PATHd8",
                                                   attempt_fix = TRUE) {
  result_matrix <- matrix(nrow = dim(patristic_matrix)[1], ncol = dim(patristic_matrix)[2])
  if (is.null(target_tree$edge.length)) {
    target_tree$edge.length <- numeric(nrow(target_tree$edge))
  }
  try(result_matrix <- phylo_to_patristic_matrix(congruify_and_check(
    reference = patristic_matrix_to_phylo(patristic_matrix),
    target = target_tree, scale = dating_method, attempt_fix = attempt_fix
  )))
  return(result_matrix)
}

# Note that originally trees were stored as patristic matrices. This was intended
# to make subsetting fast. The downside is large memory usage. Klaus Schliep wrote
# fast tree subsetting for phylo and multiphylo objects, so now trees are stored
# internally as objects of this type, but with the final output after pruning
# going through patristic matrices.

# Some trees are so large that they can't be stored as patristic distance matrices. For all others,
# patristic matrices are better. For example, for the 20,000 HeathEtAl2012 trees of 35 taxa,
# getting a subset down to two taxa takes 0.0475 seconds just for the pruning, 0.0504 seconds
# for pruning and getting a subset, for a single tree (run times go up linearly with number of trees:
# pruning and converting 1000 trees takes 3 seconds). Subsetting 1000 trees from the patristic
# distance matrix takes just 0.0013 seconds.

# in case we want to cache. Not clear we do.
# Used inside: patristic_matrix_array_phylo_congruify, phylo_get_subset_array and phylo_congruify
#' Get a patristic matrix from a `phylo` object.
#' @inheritParams phylo_check
#' @inheritParams congruify_and_check
#' @param test Default to `TRUE`. Whether to test if `phy` has branch lengths and is ultrametric or not.
#' @return A patristic matrix.
phylo_to_patristic_matrix <- function(phy, test = TRUE, tol = 0.01, option = 2) {
  # stores the distance between taxa
  patristic_matrix <- NA
  if (inherits(phy, "phylo")) {
    if (test) {
      if (!ape::is.ultrametric(phy, tol = tol, option = option)) {
        stop("Currently, datelife require that chronograms are ultrametric.") # can pad them so that terminals all reach to present
      }
    }
    patristic_matrix <- stats::cophenetic(phy)
  }
  return(patristic_matrix)
}

# Used inside: get_subset_array_dispatch.
#' Subset a reference and a target tree given as `phylo` objects.
#' @inheritParams get_subset_array_dispatch
#' @inheritParams phylo_get_subset_array
#' @inherit phylo_get_subset_array return
phylo_subset_both <- function(reference_tree,
                              taxa,
                              phy = NULL,
                              phy4 = NULL,
                              dating_method = "PATHd8") {
  # COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
  #  if (inherits(reference_tree, "phylo")) {
  #    reference_tree<-c(reference_tree) #from here in, assumes multiphylo object, even if a single tree
  #  }
  congruify <- FALSE
  if (inherits(phy, "phylo")) {
    congruify <- TRUE
  } else {
    congruify <- FALSE
  }
  if (congruify) {
    return(phylo_get_subset_array_congruify(reference_tree = reference_tree,
                                            taxa = taxa,
                                            phy = phy,
                                            dating_method = dating_method))
  } else { # when congruify is FALSE:
    return(phylo_get_subset_array(reference_tree = reference_tree,
                                  taxa = taxa,
                                  phy4 = phy4,
                                  dating_method = dating_method))
  }
}

# Used inside: phylo_subset_both, when we don't congruify
#' Get a subset array from a `phylo` object
#' @param reference_tree A `phylo` object.
#' @inheritParams get_subset_array_dispatch
#' @inherit patristic_matrix_array_subset return
phylo_get_subset_array <- function(reference_tree,
                                   taxa,
                                   phy4 = NULL,
                                   dating_method = "PATHd8") {
  final.size <- sum(reference_tree$tip.label %in% taxa) # returns number of matches
  if (final.size >= 2) { # it's worth doing the pruning
    reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa)
    # phylo_prune_missing_taxa (PruneTree before) is the new, fast fn from Klaus Schliep.
    # Eventually will be in phangorn, currently in datelife
  }
  problem <- "none"
  patristic_matrix_array <- NA
  if (final.size < length(taxa)) {
    problem <- "Missing some taxa on chronogram, so this is probably an underestimate." # fewer taxa on final matrix than we asked for
    if (final.size < 2) {
      problem <- "Insufficient species to get an MRCA (either 1 or 0)." # we either have one species or zero. Not enough for an MRCA
      patristic_matrix_array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
    patristic_matrix_array <- phylo_to_patristic_matrix(reference_tree)
  }
  if (!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "'input' of taxa are not a clade, so this is probably an overestimate."
    }
  }
  return(list(patristic_matrix_array = patristic_matrix_array, problem = problem))
}

# Used inside: phylo_subset_both.
#' Get a congruified subset array from a `phylo` object
#' @inheritParams phylo_get_subset_array
#' @inheritParams get_subset_array_dispatch
#' @inherit patristic_matrix_array_subset return
phylo_get_subset_array_congruify <- function(reference_tree,
                                             taxa,
                                             phy = NULL,
                                             dating_method = "PATHd8") {
  final.size <- sum(reference_tree$tip.label %in% taxa) # returns number of matches
  if (final.size >= 2) { # it's worth doing the pruning
    reference_tree <- phylo_prune_missing_taxa(reference_tree, taxa)
    # phylo_prune_missing_taxa (used to be names PruneTree) is the new, fast fn from Klaus Schliep.
    # Eventually will be in phangorn, currently in datelife
  }
  problem.new <- "none"
  patristic_matrix_array.new <- NA
  if (final.size < length(taxa)) {
    problem.new <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3) {
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

#' Prune missing taxa from a `phylo` object
#' Used inside phylo_get_subset_array and phylo_get_subset_array_congruify.
#' @inheritParams get_subset_array_dispatch
#' @return A `phylo` object.
phylo_prune_missing_taxa <- function(phy, taxa) {
  return(ape::drop.tip(phy, tip = phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

# Used inside: phylo_get_subset_array_congruify.
#' Congruify a reference tree and a target tree given as `phylo` objects.
#' @inheritParams summary_matrix_to_phylo
#' @inheritParams phylo_get_subset_array
#' @inheritParams congruify_and_check
#' @return A matrix.
phylo_congruify <- function(reference_tree,
                            target_tree,
                            dating_method = "PATHd8",
                            attempt_fix = TRUE) {
  result_matrix <- matrix(nrow = ape::Ntip(reference_tree), ncol = ape::Ntip(reference_tree))
  if (is.null(target_tree$edge.length)) {
    target_tree$edge.length <- numeric(nrow(target_tree$edge)) # makes it so that branches that don't match reference tree get zero length
  }
  try(result_matrix <- phylo_to_patristic_matrix(congruify_and_check(reference = reference_tree, target = target_tree, scale = dating_method, attempt_fix = attempt_fix)))
  return(result_matrix)
}

# Used inside: patristic_matrix_array_phylo_congruify and phylo_congruify.
#' Congruify and Check.
#' @inheritParams geiger::congruify.phylo
#' @inheritParams ape::is.ultrametric
#' @param attempt_fix Default to `TRUE`. If congruification results in NA branch
#' lengths, it will attempt to fix them.
congruify_and_check <- function(reference, target, taxonomy = NULL, tol = 0.01,
                                option = 2, scale = "pathd8", attempt_fix = TRUE) {
  if (!ape::is.ultrametric(reference, tol = tol, option = option)) {
    return(NA)
  }
  new.tree <- phylo_tiplabel_underscore_to_space(suppressWarnings(geiger::congruify.phylo(
    phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(target),
    taxonomy = taxonomy, tol = tol, scale = scale
  )$phy)) # suppressing warnings b/c geiger ignores tolerance
  if (anyNA(new.tree$edge.length) & attempt_fix) {
    message("Congruification resulted in NA edge lengths. Resolving polytomies and making up starting branch lengths")
    new.tree <- phylo_tiplabel_underscore_to_space(geiger::congruify.phylo(
      phylo_tiplabel_space_to_underscore(reference), phylo_tiplabel_space_to_underscore(
        ape::compute.brlen(ape::multi2di(target))
      ), taxonomy, tol, scale,
      ncores = 1
    )$phy)
    if (anyNA(new.tree$edge.length)) {
      message("There are still NAs in edge lengths; returning NA")
      new.tree <- NA
    }
  }
  new.tree$edge.length[which(new.tree$edge.length < 0)] <- 0 # sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
  return(new.tree)
}


#' Convert spaces to underscores in trees.
#'
#' `phylo_tiplabel_space_to_underscore` is used in: [make_mrbayes_runfile()],
#' [tree_get_singleton_outgroup()],
#' [congruify_and_check()], [patristic_matrix_array_phylo_congruify()].
#'
#' @inheritParams phylo_check
#' @return A `phylo` object.
phylo_tiplabel_space_to_underscore <- function(phy) {
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  return(phy)
}

#' Convert underscores to spaces in trees.
#'
#' `phylo_tiplabel_underscore_to_space` is used inside [patristic_matrix_array_phylo_congruify()], [congruify_and_check()].
#'
#' @inheritParams phylo_check
#' @return A `phylo` object.
phylo_tiplabel_underscore_to_space <- function(phy) {
  # a better name for this function would be underscore_to_blank
  # add method .phylo
  # change tip and node labels
  phy$tip.label <- gsub("_", " ", phy$tip.label)
  # make sure there is only one consecutive blank at a time
  return(phy)
}
