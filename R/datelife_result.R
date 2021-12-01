#' Get a list of patristic matrices from a given \code{datelifeQuery} object
#'
#'
#' @param input A \code{datelifeQuery} object from [make_datelife_query()].
#' @inheritParams get_datelife_result
#' @inheritDotParams make_datelife_query
#' @return A \code{datelifeResult} object - a named list of patristic matrices.
#' @details If there is just one taxon name in \code{input$cleaned_names}, the
#' function will run [make_datelife_query()] setting \code{get_spp_from_taxon = TRUE}.
#' The \code{datelifeQuery} used as input can be accessed with \code{attributes(datelifeResult)$query}.
#' @export
get_datelife_result_datelifequery <- function(input = NULL,
                                              partial = TRUE,
                                              cache = "opentree_chronograms",
                                              ...) {
  #
  message("... Getting a DateLife result.")
  if ("opentree_chronograms" %in% cache) {
    utils::data("opentree_chronograms", package = "datelife")
    cache <- get("opentree_chronograms")
  }
  if (suppressMessages(!is_datelife_query(input))) {
    stop("'input' must be a 'datelifeQuery' object.")
  }
  if (length(input$cleaned_names) == 1) {
    message("Can't get divergence times from just one taxon in 'input$cleaned_names'.")
    message("Making a DateLife Query again, with 'get_spp_from_taxon = TRUE'.")
    input <- make_datelife_query(
      input = input$cleaned_names,
      get_spp_from_taxon = TRUE,
      ...
    )
  }
  # setting phy to NULL always; it is a bad idea to congruify subset trees,
  # do that later in summarizing steps
  results_list <- lapply(cache$trees,
    get_subset_array_dispatch,
    taxa = input$cleaned_names,
    phy = NULL
  )
  datelife_result <- results_list_process(results_list,
    input$cleaned_names,
    partial = partial
  )
  message("DateLife result obtained!")
  class(datelife_result) <- c("datelifeResult")
  attr(datelife_result, "datelife_query") <- input
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

#' Get a list of patristic matrices from a given set of taxon names
#'
#' @description Go from a vector of species, newick string, or phylo object
#'
#' @inheritParams datelife_search
#' @inheritDotParams make_datelife_query
#' @export
get_datelife_result <- function(input = NULL,
                                partial = TRUE,
                                cache = "opentree_chronograms",
                                ...) {
  #
  if ("opentree_chronograms" %in% cache) {
    utils::data("opentree_chronograms", package = "datelife")
    cache <- get("opentree_chronograms")
  }
  input_dq <- input
  if (suppressMessages(!is_datelife_query(input))) {
    input_dq <- make_datelife_query(
      input = input,
      ...
    )
  }
  res <- get_datelife_result_datelifequery(
    input = input_dq,
    partial = partial,
    cache = "opentree_chronograms",
    ...
  )
}
