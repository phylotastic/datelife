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
#' 	 \item{A character vector}{With taxon names as a single comma separated
#' 	 	 starting or concatenated with [c()].}
#' 	 \item{A phylogenetic tree with taxon names as tip labels}{As a `phylo` or
#' 	 	 `multiPhylo` object, OR as a newick character string.}
#' 	\item{A `datelifeQuery` object}{An output from [make_datelife_query()].}
#' 	}
#' @inheritParams make_datelife_query
#' @inheritParams get_datelife_result_datelifequery
#' @inheritParams summarize_datelife_result
#' @details
#' If only one taxon name is given as `input`, `get_spp_from_taxon` is
#' always set to `TRUE`.
#' @inherit summarize_datelife_result return
#' @examples
#' \dontrun{ # This is a flag for package development. You are welcome to run the example.
#'
#' # For this example, we will set a temp working directory, but you can set
#' # your working directory as needed:
#' # we will use the tempdir() function to get a temporary directory:
#' tempwd <- tempdir()
#'
#' # Obtain median ages from a set of source chronograms in newick format:
#' ages <- datelife_search(c(
#'   "Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#'   "Mus musculus"
#' ), summary_format = "newick_median")
#'
#' # Save the tree in the temp working directory in newick format:
#' write(ages, file = file.path(tempwd, "some.bird.ages.txt"))
#'
#' # Obtain median ages from a set of source chronograms in phylo format
#' # Will produce same tree as above but in "phylo" format:
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
#'
#' # Save "phylo" object in newick format
#' write.tree(ages.again, file = file.path(tempwd, "some.bird.tree.again.txt"))
#'
#' # Obtain MRCA ages and target chronograms from all source chronograms
#' # Generate an htm"l output readable in any web browser:
#' ages.html <- datelife_search(c(
#'   "Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#'   "Mus musculus"
#' ), summary_format = "html")
#' write(ages.html, file = file.path(tempwd, "some.bird.trees.html"))
#' if(interactive()) system(paste("open", file.path(tempwd, "some.bird.trees.html")))
#' } # end dontrun
#' @export
datelife_search <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                            use_tnrs = FALSE,
                            # approximate_match = TRUE,
                            get_spp_from_taxon = FALSE,
                            partial = TRUE,
                            cache = "opentree_chronograms",
                            summary_format = "phylo_all",
                            na_rm = FALSE,
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
    na_rm = na_rm,
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



#' Check if we obtained an empty search with the given taxon name(s).
#' @inheritParams summarize_datelife_result
#' @inheritParams datelife_search
#' @return Boolean. If `TRUE`, no chronograms were found for the given taxon name(s).
#' If `FALSE`, the chronogram search was successful.
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
