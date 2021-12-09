#' Convert patristic matrix to a newick string. Used inside: summarize_datelife_result.
#' @param patristic_matrix A patristic matrix
#' @return A newick string
#' @export
patristic_matrix_to_newick <- function(patristic_matrix) {
  tree <- patristic_matrix_to_phylo(patristic_matrix)
  if (inherits(tree, "phylo")) {
    return(ape::write.tree(tree))
  }
  return(NA)
}
