#' Core function to generate results
#' @param input A newick string or vector of taxa
#' @param format The output format
#' @param partial How to deal with trees that have a subset of taxa in the query
#' @param plot.width Width in pixels for output plot
#' @param plot.height Height in pixels for output plot
#' @param usetnrs Whether to use OpenTree's TNRS for the input
#' @param approximatematch If using TNRS, use approximate matching
#' @param datelife.cache The list of lists containing the input trees and other info
#' @return results in the desired format
#' @export
run<-function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), format="citations", partial="yes", plot.width=600, plot.height=600, usetnrs="no", approximatematch="yes", datelife.cache=NULL) {
	if(is.null(datelife.cache)) {
	  utils::data(opentree_chronograms)
	  }
  #convert from HTML input to Boolean
  partial <- ifelse(partial=="yes", TRUE, FALSE)
  usetnrs <- ifelse(usetnrs=="yes", TRUE, FALSE)
  approximatematch <- ifelse(approximatematch=="yes", TRUE, FALSE)
  filtered.results <- GetFilteredResults(input, partial, usetnrs, approximatematch, datelife.cache)
  SummarizeResults(filtered.results, output.format=format, partial, datelife.cache, suppress.citations=TRUE)
}
