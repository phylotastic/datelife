#' Get a list of patristic matrices from a given \code{datelifeQuery} object
#'
#'
#' @param input A \code{datelifeQuery} object from
#' @return A datelifeResult object - a named list of patristic matrices. You can
#' access the original query with attributes(object)$query
#' @export
get_datelife_result_datelifequery <- function(input = NULL,
																partial = TRUE,
																cache = "opentree_chronograms"){

   if (!is_datelife_query(input)) {
   	stop("'input' must be a 'datelifeQuery' object.")
   }
	 if (length(input$cleaned_names) == 1){
		 message("Cannot get divergence times from just one taxon.")
		 if(methods::hasArg("get_spp_from_taxon")){
			 get_spp_from_taxon <- list(...)$get_spp_from_taxon
			 if(!is.logical(get_spp_from_taxon)){
				 stop("'get_spp_from_taxon' argument must be one of TRUE or FALSE")
			 }
		 } else {
			 get_spp_from_taxon <- FALSE
		 }

		input <- make_datelife_query_vector(input = input$cleaned_names,
																	 			...)
 	 }
	 if("opentree_chronograms" %in% cache){
		 utils::data("opentree_chronograms", package = "datelife")
		 cache <- get("opentree_chronograms")
	 }

   # setting phy to NULL always; it is a bad idea to congruify subset trees,
   # do that later in summarizing steps
   results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = input$cleaned_names, phy = NULL)
   datelife_result <- results_list_process(results_list, input$cleaned_names, partial)
   datelife_result_check(datelife_result, use_tnrs)
   class(datelife_result) <- c("datelifeResult")
   attr(datelife_result, "query") <- input
   return(datelife_result)
}

#' Get a list of patristic matrices from a given set of taxon names
#'
#'
#'
#' @param input A character vector of taxon names
#' @return A datelifeResult object - a named list of patristic matrices. You can
#' access the \code{datelifeQuery} object with attributes(object)$query
#' @export
get_datelife_result_vector <- function(input = NULL,
																partial = TRUE,
																cache = "opentree_chronograms",
																...){
   if (!is.character(input)) {
   	stop("'input' must be a character vector of taxon names.")
   }
	 if("opentree_chronograms" %in% cache){
		 utils::data("opentree_chronograms", package = "datelife")
		 cache <- get("opentree_chronograms")
	 }
	 if(methods::hasArg("get_spp_from_taxon")){
		 get_spp_from_taxon <- list(...)$get_spp_from_taxon
		 if(!is.logical(get_spp_from_taxon)){
			 stop("'get_spp_from_taxon' argument must be one of TRUE or FALSE")
		 }
	 } else {
		 get_spp_from_taxon <- FALSE
	 }
	 if (length(input) == 1) {
		 message("Can't get divergence times from just one taxon.")
		 if (!get_spp_from_taxon) {
			 message("Setting up get_spp_from_taxon = TRUE")
			 get_spp_from_taxon <- TRUE
		 }
	 }
	 input_dq <- make_datelife_query(input = input,
																	 get_spp_from_taxon = get_spp_from_taxon,
																	 ...)
 	 if (length(input_dq$cleaned_names) == 1) {
		 warning("'input' contains only one lineage, divergence times can't be obtained.")
		 return(NA)
 	 }

   # setting phy to NULL always; it is a bad idea to congruify subset trees,
   # do that later in summarizing steps
   results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = input$cleaned_names, phy = NULL)
   datelife_result <- results_list_process(results_list, input$cleaned_names, partial)
   datelife_result_check(datelife_result, use_tnrs)
   class(datelife_result) <- c("datelifeResult")
   attr(datelife_result, "query") <- input
   return(datelife_result)
}

#' Get a list of patristic matrices from a given set of taxon names
#'
#' @description Go from a vector of species, newick string, or phylo object
#'
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
# #' @inheritParams make_bold_otol_tree
# #' @inheritDotParams make_bold_otol_tree
#' @export
get_datelife_result <- function(input = NULL,
																partial = TRUE,
																use_tnrs = FALSE,
																approximate_match = TRUE,
																update_opentree_chronograms = FALSE,
																cache = "opentree_chronograms",
																get_spp_from_taxon = FALSE,
																verbose = FALSE){
	if(update_opentree_chronograms){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	} else {
		if("opentree_chronograms" %in% cache){
			utils::data("opentree_chronograms", package = "datelife")
			cache <- get("opentree_chronograms")
		}
	}
	if(is_datelife_query(input)){
		input_dq <- input
	} else {
		input_dq <- make_datelife_query(input = input, use_tnrs = use_tnrs,
			approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon,
			verbose = verbose)
	}

}
