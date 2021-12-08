#' Core function to generate results
#' @param input A newick string or vector of taxa
#' @param format The output format
#' @param partial How to deal with trees that have a subset of taxa in the query
#' @param plot.width Width in pixels for output plot
#' @param plot.height Height in pixels for output plot
#' @param use_tnrs Whether to use OpenTree's TNRS for the input
# # ' @param approximate_match If using TNRS, use approximate matching
#' @param opentree_chronograms The list of lists containing the input trees and other info
#' @return results in the desired format
#' @export
run <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                format = "citations",
                partial = "yes",
                plot.width = 600,
                plot.height = 600,
                use_tnrs = "no",
                # approximate_match = "yes",
                opentree_chronograms = NULL) {
  if (is.null(opentree_chronograms)) {
    utils::data(opentree_chronograms)
    opentree_chronograms <- utils::getAnywhere("opentree_chronograms")
  }
  # convert from HTML input to Boolean
  partial <- ifelse(partial == "yes", TRUE, FALSE)
  use_tnrs <- ifelse(use_tnrs == "yes", TRUE, FALSE)
  # approximate_match <- ifelse(approximate_match == "yes", TRUE, FALSE)
  datelife_result <- get_datelife_result(input,
                                         partial,
                                         use_tnrs,
                                         # approximate_match,
                                         opentree_chronograms)
  summarize_datelife_result(datelife_result = datelife_result,
    summary_format = format,
    partial,
    opentree_chronograms)
}
