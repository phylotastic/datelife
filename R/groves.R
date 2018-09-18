#' Function for computing n overlap for two vectors of names (ie., phy1$tip.label, phy2$tip.label) and seeing if they have n overlap
#'
#' n-overlap comes from Definition 2.8 of Ané et al. 10.1007/s00026-009-0017-x Groves of Phylogenetic Trees
#' @param names_1 First vector of names
#' @param names_2 Second vector of names
#' @param n Degree of overlap required
#' @return Boolean for whether the degree of overlap was met
is_n_overlap <- function(names_1, names_2, n=2) {
  return(sum(names_1 %in% names_2) >= n)
}

#' Build grove matrix
#'
#' Using theorem 1.1 of Ané et al. 10.1007/s00026-009-0017-x Groves of Phylogenetic Trees
#' @param datelife_result datelifeResult object (named list of patristic matrices)
#' @param n Degree of overlap required
#' @return matrix; each cell shows whether n-overlap exists between that pair of inputs
build_groves_matrix <- function(datelife_result, n=2) {
  grove_matrix <- matrix(FALSE, nrow=length(datelife_result), ncol=length(datelife_result))
  for (i in sequence(length(datelife_result))) {
    for (j in i:length(datelife_result)) {
        grove_matrix[i,j] <- is_n_overlap(names_1 = rownames(datelife_result[[i]]),
        names_2 = rownames(datelife_result[[j]]), n = n)
    }
  }
  return(grove_matrix)
}

#' Build grove list
#'
#' Using theorem 1.1 of Ané et al. 10.1007/s00026-009-0017-x Groves of Phylogenetic Trees
#' @param datelife_result datelifeResult object (named list of patristic matrices)
#' @param n Degree of overlap required
#' @return list of vectors; each list element is a grove
build_grove_list <- function(datelife_result, n = 2) {
  grove_matrix <- build_groves_matrix(datelife_result, n)
  grove_list <- list()
  for (i in sequence(nrow(grove_matrix))) {
    if (i == 1) {
      grove_list <- list(c(which(grove_matrix[i,])))
    } else {
      elements <- which(grove_matrix[i,])
      matching_grove <- NULL
      for (grove_index in sequence(length(grove_list))) {
        if(any(elements %in% grove_list[[grove_index]])) {
          matching_grove <- grove_index
          break()
        }
      }
      if(!is.null(matching_grove)) {
        grove_list[[matching_grove]] <- unique(c(grove_list[[matching_grove]], elements))
      } else {
        grove_list[[length(grove_list)+1]] <- elements
      }
    }
  }
  return(grove_list)
}

#' Find the grove in the case of multiple groves
#'
#' @param grove_list A list of vectors of tree indices. Each element is a grove
#' @param criterion Whether to get the grove with the most trees or the most taxa
#' @param datelife_result datelifeResult object (named list of patristic matrices). Only needed for "taxa" criterion
#' @return a vector of the elements of the chosen grove
pick_grove <- function(grove_list, criterion=c("trees", "taxa"), datelife_result) {
  criterion <- match.arg(criterion, c("trees", "taxa"))
  if(length(grove_list)==1) {
    return(grove_list[[1]])
  } else {
    if(criterion=="trees") {
      tree_counts <- lapply(grove_list, length)
      return(grove_list[[which.max(tree_counts)]])
    } else {
      taxa_counts <- rep(0, length(grove_list))
      for (grove_index in sequence(length(grove_list))) {
        taxa_in_grove <- c()
        pointers_to_trees <- grove_list[[grove_index]]
        for (tree_index in sequence(length(grove_list[[grove_index]]))) {
          taxa_in_grove <- unique(c(taxa_in_grove, rownames(datelife_result[[ pointers_to_trees[tree_index] ]])))
        }
        taxa_counts[grove_index] <- length(taxa_in_grove)
      }
      return(grove_list[[which.max(taxa_counts)]])
    }
    stop("If you are here, something went wrong in pick_grove")
  }
}

#' Filter a datelifeResult to find the largest grove
#' @param datelife_result datelifeResult object (named list of patristic matrices). Only needed for "taxa" criterion
#' @param criterion Whether to get the grove with the most trees or the most taxa
#' @param n Degree of overlap required
#' @return A datelifResult object filtered to only include one grove of trees
filter_for_grove <- function(datelife_result, criterion=c("trees", "taxa"), n = 2) {
  grove_list <- build_grove_list(datelife_result, n)
  final_trees <- pick_grove(grove_list, criterion, datelife_result)
  return(datelife_result[final_trees])
}
