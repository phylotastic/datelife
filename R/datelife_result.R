#' Get a list of patristic matrices from a given `datelifeQuery` object
#'
#' @inheritParams get_taxon_summary
#' @param partial Whether to return or exclude partially matching source chronograms,
#'   i.e, those that match some and not all of taxa given in `datelife_query`.
#'   Options are `TRUE` or `FALSE`. Defaults to `TRUE`: return all matching source
#'   chronograms.
#' @param cache A character vector of length one, with the name of the data object
#'   to cache. Default to `"opentree_chronograms"`, a data object storing Open Tree of
#'   Life's database chronograms and other associated information.
#' @param update_opentree_chronograms Whether to update the chronogram database or not.
#'   Defaults to `FALSE`.
#' @inheritDotParams make_datelife_query
#' @return A `datelifeResult` object -- a named list of patristic matrices.
#' @details If there is just one taxon name in `input$cleaned_names`, the
#' function will run [make_datelife_query()] setting `get_spp_from_taxon = TRUE`.
#' The `datelifeQuery` used as input can be accessed with `attributes(datelifeResult)$query`.
#' @export
get_datelife_result_datelifequery <- function(datelife_query = NULL,
                                              partial = TRUE,
                                              cache = "opentree_chronograms",
                                              update_opentree_chronograms = FALSE,
                                              ...) {
  #
  if (update_opentree_chronograms) {
    cache <- update_datelife_cache(write = TRUE)
  } else {
    if ("opentree_chronograms" %in% cache) {
      utils::data("opentree_chronograms", package = "datelife")
      cache <- get("opentree_chronograms")
    }
  }
  message("... Searching DateLife's OpenTree chronogram database version v",
          cache$version)
  if (suppressMessages(!is_datelife_query(input = datelife_query))) {
    stop("'datelife_query' must be a 'datelifeQuery' object.")
  }
  if (length(datelife_query$cleaned_names) == 1) {
    message("Can't get divergence times from just one taxon in 'datelife_query$cleaned_names'.")
    message("Making a 'datelifeQuery' again, setting 'get_spp_from_taxon = TRUE'.")
    datelife_query <- make_datelife_query(input = datelife_query$cleaned_names,
                                          get_spp_from_taxon = TRUE,
                                          ...)
  }
  # setting phy to NULL always; it is a bad idea to congruify subset trees,
  # do that later in summarizing steps
  results_list <- lapply(cache$trees,
    get_subset_array_dispatch,
    taxa = datelife_query$cleaned_names,
    phy = NULL
  )
  datelife_result <- results_list_process(results_list,
    datelife_query$cleaned_names,
    partial = partial
  )
  message("Search done!")
  message("\nInput taxon names were found in ", length(datelife_result), " chronograms.")
  class(datelife_result) <- c("datelifeResult")
  attr(datelife_result, "datelife_query") <- datelife_query
  return(datelife_result)
}

# The following function is a get_datelife_result method that takes vectors as input
# It might not be necessary
# #' Get a list of patristic matrices from a given set of taxon names
# #'
# #'
# #'
# #' @param input A character vector of taxon names
# #' @inheritParams get_datelife_result_datelifequery
# #' @return A datelifeResult object - a named list of patristic matrices. You can
# #' access the \code{datelifeQuery} object with attributes(object)$query
# #' @export
# get_datelife_result_vector <- function(input = NULL,
# 																partial = TRUE,
# 																cache = "opentree_chronograms",
# 																...){
#    if (!is.character(input)) {
#    	stop("'input' must be a character vector of taxon names.")
#    }
# 	 if("opentree_chronograms" %in% cache){
# 		 utils::data("opentree_chronograms", package = "datelife")
# 		 cache <- get("opentree_chronograms")
# 	 }
# 	 if(methods::hasArg("get_spp_from_taxon")){
# 		 get_spp_from_taxon <- list(...)$get_spp_from_taxon
# 		 if(!is.logical(get_spp_from_taxon)){
# 			 stop("'get_spp_from_taxon' argument must be one of TRUE or FALSE")
# 		 }
# 	 } else {
# 		 get_spp_from_taxon <- FALSE
# 	 }
# 	 if (length(input) == 1) {
# 		 message("Can't get divergence times from just one taxon.")
# 		 if (!get_spp_from_taxon) {
# 			 message("Setting up get_spp_from_taxon = TRUE")
# 			 get_spp_from_taxon <- TRUE
# 		 }
# 	 }
# 	 input_dq <- make_datelife_query(input = input,
# 																	 get_spp_from_taxon = get_spp_from_taxon,
# 																	 ...)
#  	 if (length(input_dq$cleaned_names) == 1) {
# 		 warning("'input' contains only one lineage, divergence times can't be obtained.")
# 		 return(NA)
#  	 }
#
#    # setting phy to NULL always; it is a bad idea to congruify subset trees,
#    # do that later in summarizing steps
#    results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = input$cleaned_names, phy = NULL)
#    datelife_result <- results_list_process(results_list, input$cleaned_names, partial)
#    class(datelife_result) <- c("datelifeResult")
#    attr(datelife_result, "datelife_query") <- input
#    return(datelife_result)
# }

#' Get a patristic matrix of time of lineage divergence data for a given set of taxon names
#'
#' `get_datelife_result` takes as input a vector of taxon names, a newick string,
#' a `phylo` object, or a`datelifeQuery` object. It searches the chronogram
#' database specified in `cache` for chronograms matching two or more given
#' taxon names. For each matching chronogram, it extracts time of lineage
#' divergence data and stores it as a patristic matrix. It then lists all
#' resulting patristic matrices. Each list element is named with the study
#' citation of the source chronogram.
#'
#' @inheritParams datelife_search
# datelife_search has param input
#' @inheritParams get_datelife_result_datelifequery
#' @inheritDotParams make_datelife_query
#' @inherit get_datelife_result_datelifequery return
#' @export
get_datelife_result <- function(input = NULL,
                                partial = TRUE,
                                cache = "opentree_chronograms",
                                update_opentree_chronograms = FALSE,
                                ...) {
  #
  input_dq <- input
  if (suppressMessages(!is_datelife_query(input))) {
    input_dq <- make_datelife_query(input = input,
                                    ...)
  }
  res <- get_datelife_result_datelifequery(datelife_query = input_dq,
                    partial = partial,
                    cache = cache,
                    update_opentree_chronograms = update_opentree_chronograms)
}
