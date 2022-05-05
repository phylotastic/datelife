#' Function for computing n-overlap for two vectors of names (ie., phy1$tip.label, phy2$tip.label) and seeing if they have n overlap
#'
#' @description This function implements definition 2.8 for n-overlap from Ané
#'   et al. (2009) \doi{10.1007/s00026-009-0017-x}.
#' @param names_1 First vector of names
#' @param names_2 Second vector of names
#' @param n Degree of overlap required
#' @return Boolean for whether the degree of overlap was met or not.
#' @references
#' Ané, C., Eulenstein, O., Piaggio-Talice, R., & Sanderson, M. J. (2009).
#' "Groves of phylogenetic trees". Annals of Combinatorics, 13(2), 139-167,
#' \doi{10.1007/s00026-009-0017-x}.
#' @export
is_n_overlap <- function(names_1, names_2, n = 2) {
  return(sum(names_1 %in% names_2) >= n)
}

#' Find the grove for a group of chronograms and build a matrix.
#'
#' @description This function implements theorem 1.1 of Ané et al. (2009) \doi{10.1007/s00026-009-0017-x}
#'   to find a grove for a given group of chronograms.
#' @param datelife_result A `datelifeResult` object.
#' @param n The degree of taxon name overlap among input chronograms. Defaults
#'   to `n = 2`, i.e., at least two overlapping taxon names.
#' @return A matrix. Each cell shows whether n-overlap exists between a pair of inputs.
#' @references
#' Ané, C., Eulenstein, O., Piaggio-Talice, R., & Sanderson, M. J. (2009).
#' "Groves of phylogenetic trees". Annals of Combinatorics, 13(2), 139-167,
#' \doi{10.1007/s00026-009-0017-x}.
#' @export
build_grove_matrix <- function(datelife_result, n = 2) {
  grove_matrix <- matrix(FALSE, nrow = length(datelife_result), ncol = length(datelife_result))
  for (i in sequence(length(datelife_result))) {
    for (j in i:length(datelife_result)) {
      grove_matrix[i, j] <- is_n_overlap(
        names_1 = rownames(datelife_result[[i]]),
        names_2 = rownames(datelife_result[[j]]), n = n
      )
    }
  }
  return(grove_matrix)
}

#' Build grove list
#'
#' @inherit build_grove_matrix description
#' @inheritParams build_grove_matrix
#' @return A list of vectors; each list element is a grove.
#' @export
build_grove_list <- function(datelife_result, n = 2) {
  grove_matrix <- build_grove_matrix(datelife_result, n)
  grove_list <- list()
  for (i in sequence(nrow(grove_matrix))) {
    if (i == 1) {
      grove_list <- list(c(which(grove_matrix[i, ])))
    } else {
      elements <- which(grove_matrix[i, ])
      matching_grove <- NULL
      for (grove_index in sequence(length(grove_list))) {
        if (any(elements %in% grove_list[[grove_index]])) {
          matching_grove <- grove_index
          break()
        }
      }
      if (!is.null(matching_grove)) {
        grove_list[[matching_grove]] <- unique(c(grove_list[[matching_grove]], elements))
      } else {
        grove_list[[length(grove_list) + 1]] <- elements
      }
    }
  }
  return(grove_list)
}

#' Pick a grove in the case of multiple groves in a set of trees.
#'
#' @param grove_list A list of vectors of tree indices. Each element is a grove.
#' @inheritParams summarize_datelife_result
#' @param datelife_result A `datelifeResult` object. Only needed for `criterion = "taxa"`.
#' @return A numeric vector of the elements of the picked grove.
#' @export
pick_grove <- function(grove_list, criterion = "taxa", datelife_result) {
  criterion <- match.arg(criterion, c("trees", "taxa"))
  if (length(grove_list) == 1) {
    return(grove_list[[1]])
  } else {
    if (criterion == "trees") {
      tree_counts <- lapply(grove_list, length)
      return(grove_list[[which.max(tree_counts)]])
    } else {
      taxa_counts <- rep(0, length(grove_list))
      for (grove_index in sequence(length(grove_list))) {
        taxa_in_grove <- c()
        pointers_to_trees <- grove_list[[grove_index]]
        for (tree_index in sequence(length(grove_list[[grove_index]]))) {
          taxa_in_grove <- unique(c(taxa_in_grove, rownames(datelife_result[[pointers_to_trees[tree_index]]])))
        }
        taxa_counts[grove_index] <- length(taxa_in_grove)
      }
      return(grove_list[[which.max(taxa_counts)]])
    }
    stop("If you are here, something went wrong in pick_grove")
  }
}

#' Filter a `datelifeResult` object to find the largest grove.
#' @inheritParams pick_grove
#' @inheritParams build_grove_matrix
#' @return A `datelifeResult` object filtered to only include one grove of trees.
#' @export
filter_for_grove <- function(datelife_result, criterion = "taxa", n = 2) {
  criterion <- match.arg(criterion, c("trees", "taxa"))
  grove_list <- build_grove_list(datelife_result, n)
  final_trees <- pick_grove(grove_list, criterion, datelife_result)
  return(datelife_result[final_trees])
}

#' Get grove from a \code{datelifeResult} object that can be converted to phylo
#' from a median summary matrix
#' @inheritParams filter_for_grove
#' @return A list of two elements:
#' \describe{
#' 	\item{best_grove}{A `datelifeResult` object filtered to only include one grove of trees that can be summarized with median or sdm.
#' 	}
#' 	\item{overlap}{The degree of taxon names overlap among trees in the best grove.
#' 	}
#' }
#' @export
# #' @details
get_best_grove <- function(datelife_result, criterion = "taxa", n = 2) {
  # for testing:
  # utils::data(subset2_taxa)
  # spp_query <- make_datelife_query(subset2_taxa)
  # datelife_result <- get_datelife_result(spp_query)
  # datelife_result <- check_datelife_result(datelife_result)
  median_nj <- NULL
  while (!inherits(median_nj, "phylo")) {
    message(paste0("Trying with overlap = ", n, "\n"))
    best_grove <- filter_for_grove(datelife_result, criterion = criterion, n = n)
    # length(best_grove)
    # we use patristic_matrix_to_phylo as a test that the grove can be clustered into a tree
    # for that we first get the median matrix and then try to cluster (with njs) catching the error
    # until we get a tree
    median_matrix <- datelife_result_median_matrix(best_grove)
    median_nj <- tryCatch(suppressMessages(suppressWarnings(ape::njs(median_matrix))),
      error = function(e) NA
    )
    # issue: sometimes max(ape::branching.times) is off (too big or too small), so we could
    # standardize by real median of original data (max(mrcas)).
    # median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
    # We might have solved the above issue by using our method developped for sdm matrices
    n <- n + 1
  }
  message(paste0("Success!", "\n"))
  class(best_grove) <- class(datelife_result)
  return(list(best_grove = best_grove, overlap = n - 1))
}
