#' Get a median summary chronogram from a `datelifeResult` object.
#' @inheritParams summarize_datelife_result
#' @inheritDotParams summary_matrix_to_phylo
#' @return A `phylo` object.
#' @export
datelife_result_median <- function(datelife_result, ...) {
  # for debugging here:
  # datelife_result <- get_best_grove(subset2_result)$best_grove
  # datelife_result <- best_grove
  message("... Calculating a median summary chronogram.")
  median_matrix <- datelife_result_median_matrix(datelife_result)
  # tree <- suppressWarnings(suppressMessages(patristic_matrix_to_phylo(median.matrix,
  # 		  clustering_method = "nj", fix_negative_brlen = TRUE)))
  # sometimes max(ape::branching.times) is off (too big or too small), so we could
  # standardize by real median of original data (max(mrcas)).
  # median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
  phy <- summary_matrix_to_phylo(median_matrix, use = "median", ...)
  phy$data <- datelife_result
  phy$citation <- names(datelife_result)
  return(phy)
}

#' Compute a median matrix of a `datelifeResult` object.
#' @inheritParams summarize_datelife_result
#' @return A patristic distance summary matrix from a `datelifeResult` object.
#' @export
datelife_result_median_matrix <- function(datelife_result) {
  # datelife_result <- check_datelife_result(datelife_result)
  patristic.array <- patristic_matrix_list_to_array(datelife_result)
  median.matrix <- summary_patristic_matrix_array(patristic.array)
  # when matrix comes from median, upgma gives much older ages than expected
  # we used nj to cluster in this case
  # now we prefer our algorithm
  return(median.matrix)
}

#' Compute a variance matrix of a `datelifeResult` object.
#' @inheritParams summarize_datelife_result
#' @return A variance matrix from a `datelifeResult` object.
#' @export
datelife_result_variance_matrix <- function(datelife_result) {
  # datelife_result <- check_datelife_result(datelife_result)
  patristic.array <- patristic_matrix_list_to_array(datelife_result)
  var.matrix <- summary_patristic_matrix_array(patristic.array, fn = stats::var)
  return(var.matrix)
}
