#' Reconstruct a supertree from a `datelifeResult` object using the Super Distance Matrix (SDM) method.
#'
#'
#'
#' @inheritParams summarize_datelife_result
#' @inheritParams make_sdm
#' @inheritDotParams summary_matrix_to_phylo
#' @return A supertree with branch lengths proportional to time, obtained by
#' summarizing individual chronograms given as input in `datelife_result`.
#' It is returned as an object of class `datelifeSDM`, which is a `phylo` object
#' with an additional `$data` element storing the input chronograms as a
#' `datelifeResult` object, and a `$citation` element containing
#' citations of studies from input chronograms.
#' @details
#' Chronograms given as input in `datelife_result` are summarized with the Super Distance
#' Matrix (SDM) method described in Criscuolo et al. (2006) \doi{10.1080/10635150600969872},
#' implemented with the function [ape::SDM()]. The resulting summary SDM is
#' clustered with [summary_matrix_to_phylo()].
#' @export
#' @references
#' Criscuolo A, Berry V, Douzery EJ, Gascuel O.
#' (2006) "SDM: a fast distance-based approach for (super) tree building in
#' phylogenomics" \doi{10.1080/10635150600969872}.

datelife_result_sdm_phylo <- function(datelife_result,
                                weighting = "flat",
                                ...) {
  # add a check for datelife_result?
  phy <- NA
  # used.studies <- names(datelife_result)
  if (length(datelife_result) == 1) {
    phy <- patristic_matrix_to_phylo(datelife_result, fix_negative_brlen = FALSE)
    unpadded.matrices <- datelife_result
  } else {
    # datelife_result <- best_grove
    unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
    good.matrix.indices <- get_goodmatrices(unpadded.matrices)
    if (length(good.matrix.indices) > 1) {
      unpadded.matrices <- unpadded.matrices[good.matrix.indices]
      SDM.result <- make_sdm(unpadded.matrices, weighting)
      # it is important to use upgma as clustering method; nj produces much younger ages when the matrix comes from sdm
      # phy <- patristic_matrix_to_phylo(SDM.result, clustering_method = clustering_method, fix_negative_brlen = TRUE)
      phy <- summary_matrix_to_phylo(SDM.result, ...) # this also contains the age distributions and calibrations used
    } else {
      if (length(good.matrix.indices) == length(datelife_result)) {
        warning("There are not enough input chronograms to run SDM. You can help by uploading trees to Open Tree of Life tree store.")
      } else {
        warning("All input chronograms throw an error when running SDM. This is not your fault.")
      }
      stop("SDM tree cannot be generated with this set of taxa.")
    }
    class(unpadded.matrices) <- "datelifeResult"
  }
  phy$data <- unpadded.matrices
  phy$citation <- names(unpadded.matrices)
  class(phy) <- c(class(phy), "datelifeSDM")
  return(phy)
}
#' Function to remove missing taxa from a `datelifeResult` object.
#' @description Used in [datelife_result_sdm_phylo()].
#' @param patristic_matrix A patristic matrix with row and column names for taxa
#' @return patristic_matrix for all_taxa
#' @export
patristic_matrix_unpad <- function(patristic_matrix) {
  bad.ones <- which(apply(is.na(patristic_matrix), 2, all))
  if (length(bad.ones) > 0) {
    patristic_matrix <- patristic_matrix[-bad.ones, -bad.ones]
  }
  return(patristic_matrix)
}

#' Make a Super Distance Matrix (SDM) from a list of good matrices obtained with [get_goodmatrices()]
#' @param unpadded.matrices A list of patristic matrices, a `datelifeResult` object.
#' @param weighting A character vector indicating how much weight to give to each
#' tree in `input` during the SDM analysis. Options are:
#' \describe{
#' 	\item{weighting = "flat"}{All trees have equal weighting.}
#' 	\item{weighting = "taxa"}{Weight is proportional to number of taxa.}
#' 	\item{weighting = "inverse"}{Weight is proportional to 1 / number of taxa.}
#' }
#' Defaults to `weighting = "flat"`.
#' @return A matrix.
#' @export
make_sdm <- function(unpadded.matrices, weighting = "flat") {
  # used.studies <- used.studies[good.matrix.indices]
  weights <- rep(1, length(unpadded.matrices))
  if (weighting == "taxa") {
    weights <- unname(sapply(unpadded.matrices, dim)[1, ])
  }
  if (weighting == "inverse") {
    weights <- 1 / unname(sapply(unpadded.matrices, dim)[1, ])
  }
  message(cat("\n", "Synthesizing", length(unpadded.matrices), "chronograms with SDM"))
  SDM.result <- do.call(ape::SDM, c(unpadded.matrices, weights))[[1]]
  return(SDM.result)
}
#' Get indices of good matrices to apply Super Distance Matrix (SDM) method with [make_sdm()].
#' @inheritParams make_sdm
#' @return A numeric vector of good matrix indices in unpadded.matrices.
#' @export
get_goodmatrices <- function(unpadded.matrices) {
  good.matrix.indices <- c()
  for (i in sequence(length(unpadded.matrices))) {
    test.result <- NA
    # Rationale here: some chronograms always cause errors with SDM, even when trying to get a consensus of them
    # with themselves. For now, throw out of synthesis.
    try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent = TRUE)
    message(cat(i, "out of", length(unpadded.matrices), "chronograms tried: "), appendLF = FALSE)
    if (is.finite(test.result)) {
      good.matrix.indices <- append(good.matrix.indices, i)
      message(cat(" Ok."))
    } else {
      message(cat(" Failed."))
    }
  }
  return(good.matrix.indices)
}
#' Go from a `datelifeResult` object to a Super Distance Matrix (SDM) using weighting = "flat"
#' @inheritParams summarize_datelife_result
#' @return A numeric matrix.
#' @export
datelife_result_sdm_matrix <- function(datelife_result) {
  unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices)
  if (length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- make_sdm(unpadded.matrices, weighting = "flat")
  }
  return(sdm_matrix)
}
