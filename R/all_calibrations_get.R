#' Get secondary calibrations available for a set of taxa.
#'
#' \code{get_all_calibrations} returns secondary calibrations for each pair of given
#' taxon names.
#'
#' @details If input is a character vector of taxon names, the function calls a
#' \code{\link[=datelife_search]{datelife_search}} with \code{summary_format}
#' argument set to \code{"phylo_all"}, to get all chronograms containing at
#' least 2 of the taxa given in \code{input}.
#' The function retrieves divergence times for each pair of taxa in any given tree with
#' branch lengths proportional to time (chronogram).
#' @param input A character vector of taxon names; a newick string, OR a phylo OR
#' multiPhylo object with branch lengths proportional to time; OR a datelifeResult object.
#' @inheritDotParams datelife_search
#' @param each Boolean, default to \code{FALSE}: all calibrations are returned in
#' the same data frame. If \code{TRUE}, calibrations from each chronogram are returned
#' in separate data frames.
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
get_all_calibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                                 ...,
																 each = FALSE) {
  # TODO: is_datelife_search_input function or any type of input format checking function
  # to replace the following long conditional
  # and trap the case were input is a list
	if (!inherits(input, "datelifeResult") & !inherits(input, "phylo") & !inherits(input, "multiPhylo")) {
		datelife_phylo <- datelife_search(input = input,
																			summary_format = "phylo_all",
																			...)
	}
	# inherits(datelife_phylo, "datelifeResult")
	if (inherits(input, "datelifeResult")) {
		datelife_phylo <- suppressMessages(summarize_datelife_result(datelife_result = input,
		                                                             summary_format = "phylo_all"))
	}
	if (inherits(input, "multiPhylo")) {
		datelife_phylo <- input
		xx <- sapply(datelife_phylo, "[", "edge.length")
		xx <- unname(sapply(xx, is.null))
	  if (all(xx)) {
	    message("trees in 'multiPhylo' input have no branch lengths. \n There are no calibrations to return.")
	    return(NA)
	  }
		if (any(xx)) {
	    ii <- which(xx)
	    message("Some trees in 'multiPhylo' input have no branch lengths.")
	    message("Will leave tree ",
	            paste(ii, collapse = " - "),
	            " out of the analysis")
	    datelife_phylo <- datelife_phylo[which(!xx)]
	  }
	}
	if (inherits(input, "phylo")) {
	  if (is.null(input$edge.length)) {
	    stop("input tree has no branch lengths")
	  }
		datelife_phylo <- list(input)
	}

	if (each) {
		# we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
	  constraints_df <- vector(mode = "list")
	} else {
	  # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
		constraints_df <- data.frame()
	}
	for (i in seq(length(datelife_phylo))) {
		local_df <- suppressWarnings(
			          geiger::congruify.phylo(reference = datelife_phylo[[i]],
				                                target = datelife_phylo[[i]],
															          scale = NA,
															          ncores = 1))$calibrations
		# suppressedWarnings bc of warning message when running
		# geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA)
		# 		Warning message:
		# In if (class(stock) == "phylo") { :
		# the condition has length > 1 and only the first element will be used
		local_df$reference <- names(datelife_phylo)[i]
		if (each) {
			constraints_df <- c(constraints_df, list(local_df))
			} else {
				if (i == 1) {
					constraints_df <- local_df
				} else {
					constraints_df <- rbind(constraints_df, local_df)
				}
		 }
	}
	if (each) {
		names(constraints_df) <- names(datelife_phylo)
	}
	attr(constraints_df, "chronograms") <- datelife_phylo
	return(constraints_df)
}
