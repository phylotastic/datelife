#' Get scientific, peer-reviewed information on time of lineage
#' divergence openly available for a given set of taxon names
#'
#' @description `datelife_search` is the core DateLife function to find and
#'   get all openly available, peer-reviewed scientific information on time of
#'   lineage divergence for a set of `input` taxon names given as a character
#'   vector, a newick character string, a `phylo` or `multiPhylo` object or as a
#'   an already processed `datelifeQuery` object obtained with [make_datelife_query()].
#'
#' @aliases datelife
#' @param input One of the following:
#' \describe{
#' 	 \item{A character vector}{With taxon names as a single comma separated strting or concatenated with [c()].}
#' 	 \item{A phylogenetic tree with taxon names as tip labels}{As a `phylo` or `multiPhylo`
#' 	 			object, OR as a newick character string.}
#' 	\item{A `datelifeQuery` object}{An output from [make_datelife_query()].}
#' 	}
#' @inheritParams make_datelife_query
#' @inheritParams get_datelife_result_datelifequery
#' @inheritParams summarize_datelife_result
#' @export
#' @details
#' If only one taxon name is given as `input`, `get_spp_from_taxon` is
#' always set to `TRUE`.
#'
#' @examples
#' # obtain median ages from a set of source chronograms in newick format:
#' ages <- datelife_search(c(
#'   "Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#'   "Mus musculus"
#' ), summary_format = "newick_median")
#' # save the tree in newick format
#' write(ages, file = "some.bird.ages.txt")
#'
#' # obtain median ages from a set of source chronograms in phylo format
#' # will produce same tree as above but in r phylo format:
#' ages.again <- datelife_search(c(
#'   "Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#'   "Mus musculus"
#' ), summary_format = "phylo_median")
#' plot(ages.again)
#' library(ape)
#' ape::axisPhylo()
#' mtext("Time (million years ago)", side = 1, line = 2, at = (max(get("last_plot.phylo",
#'   envir = .PlotPhyloEnv
#' )$xx) * 0.5))
#' write.tree(ages.again, file = "some.bird.tree.again.txt") # saves phylo object in newick format
#'
#' # obtain mrca ages and target chronograms from all source chronograms
#' # generate an html  output readable in any web browser:
#' ages.html <- datelife_search(c(
#'   "Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#'   "Mus musculus"
#' ), summary_format = "html")
#' write(ages.html, file = "some.bird.trees.html")
#' system("open some.bird.trees.html")
datelife_search <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                            use_tnrs = FALSE,
                            # approximate_match = TRUE,
                            get_spp_from_taxon = FALSE,
                            partial = TRUE,
                            cache = "opentree_chronograms",
                            summary_format = "phylo_all",
                            na.rm = FALSE,
                            summary_print = c("citations", "taxa"),
                            taxon_summary = c("none", "summary", "matrix"),
                            criterion = "taxa") {
  # TODO: find a way not to repeat partial and cache arguments, they are used in
  # datelife_search,  get_datelife_result and summarize_datelife_result
  # maybe take out option to update cache and work with versions of cache for reproducibility
  message("... Running a DateLife search.")
  datelife_query <- input
  if (suppressMessages(!is_datelife_query(input))) {
    datelife_query <- make_datelife_query(
      input = input,
      use_tnrs = use_tnrs,
      # approximate_match = approximate_match,
      get_spp_from_taxon = get_spp_from_taxon
    )
  }
  datelife_result.here <- get_datelife_result_datelifequery(
    datelife_query = datelife_query,
    partial = partial,
    # approximate_match = approximate_match,
    cache = cache)
  # print.datelife(datelife_result = datelife_result.here)
  # datelife <- list(datelife_query = datelife_query, datelife_result = datelife_result.here)
  # class(datelife) <- "datelife"
  # return(datelife)
  res <- summarize_datelife_result(
    datelife_result = datelife_result.here,
    datelife_query = datelife_query,
    summary_format = summary_format,
    na.rm = na.rm,
    summary_print = summary_print,
    taxon_summary = taxon_summary,
    criterion = criterion
  )
  attr(res, "datelife_result") <- datelife_result.here
  message("DateLife search done!")
  return(res)
}

# print.datelife <- function(datelife){
# 	cat("Number and which queried taxa found in chronogram database")
# 	cat("Number of queried taxa not found")
# 	cat("Number of chronograms with at least two queried taxa found in database")
# }



#' Check if we obtained an empty search with the set of input taxon names
#' @inheritParams summarize_datelife_result
#' @inheritParams datelife_search
#' @export
is_datelife_result_empty <- function(datelife_result,
                                     use_tnrs = FALSE) {
  if (length(datelife_result) < 1) {
    warning("'datelife_result' object is empty.", call. = FALSE)
    message("'input' taxon names were not found in the local chronogram database.")
    if (!use_tnrs) {
      message("setting 'use_tnrs = TRUE' might change this, but it is time consuming.")
    }
    return(TRUE)
  }
  return(FALSE)
}

#' \code{datelifeResult} object of three birds "Rhea americana", "Pterocnemia pennata", and "Struthio camelus"
#'
#' @name threebirds_dr
#' @docType data
#' @format A list of 9 named patristic matrix
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree subset chronogram ants datelife
#' @details
#' Generated with:
#' threebirds_dr <- get_datelife_result(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
#' partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, cache = "opentree_chronograms")
#' use_data(threebirds_dr)
"threebirds_dr"
