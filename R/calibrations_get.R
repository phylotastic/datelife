#' Extract secondary calibrations from a \code{phylo} or \code{multiPhylo}
#' object with branch lengths proportional to time
#'
#' @description \code{extract_calibrations_phylo} extracts divergence times
#' (i.e., secondary calibrations) for each taxon pair in a given
#' \code{phylo} or \code{multiPhylo} object.
#'
#' @param input a \code{phylo} or \code{multiPhylo} object with branch lengths
#' proportional to time.
#' @param each Boolean, default to \code{FALSE}: all calibrations are returned in
#' the same data frame. If \code{TRUE}, calibrations from each chronogram are returned
#' in separate data frames.
#' @return An object of class \code{datelifeCalibrations} -- a \code{data frame}
#' of secondary calibrations (or list of \code{data frames}, if \code{each = TRUE}),
#' for each pair of taxon names in \{input]. The attribute "chronograms" stores
#' the source data from which the calibrations were extracted.
#' @export
extract_calibrations_phylo <- function(input = NULL,
																       each = FALSE) {
	if (inherits(input, "multiPhylo")) {
		chronograms <- input
		xx <- sapply(chronograms, "[", "edge.length")
		xx <- unname(sapply(xx, is.null))
	  if (all(xx)) {
	    warning("trees in 'multiPhylo' input have no branch lengths. \n There are no calibrations to return!")
	    return(NA)
	  }
		if (any(xx)) {
	    ii <- which(xx)
	    message("Some trees in 'multiPhylo' input have no branch lengths.")
	    message("Will leave tree(s) ",
	            paste(ii, collapse = " - "),
	            " out of the analysis.")
	    chronograms <- chronograms[which(!xx)]
	  }
	}
	if (inherits(input, "phylo")) {
	  if (is.null(input$edge.length)) {
			warning("'input' tree has no branch lengths. \n There are no calibrations to return!")
	    return(NA)
	  }
		chronograms <- list(input)
	}

	if (each) {
	  calibrations <- vector(mode = "list")
	} else {
	  # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
		calibrations <- data.frame()
	}
	for (i in seq(length(chronograms))) {
		chronograms[[i]]$tip.label <- gsub(" ", "_", chronograms[[i]]$tip.label) # the battle will never end!
		local_df <- suppressWarnings(
			          geiger::congruify.phylo(reference = chronograms[[i]],
				                                target = chronograms[[i]],
															          scale = NA,
															          ncores = 1))$calibrations
		# suppressedWarnings bc of warning message when running
		# geiger::congruify.phylo(reference = chronograms[[i]], target = chronograms[[i]], scale = NA)
		# 		Warning message:
		# In if (class(stock) == "phylo") { :
		# the condition has length > 1 and only the first element will be used
		local_df$reference <- names(chronograms)[i]
		if (each) {
			calibrations <- c(calibrations, list(local_df))
			} else {
				if (i == 1) {
					calibrations <- local_df
				} else {
					calibrations <- rbind(calibrations, local_df)
				}
		 }
	}
	if (each) {
		names(calibrations) <- names(chronograms)
	}
	attr(calibrations, "chronograms") <- chronograms
	# TODO check that class data frame is also preserved. Might wanna do:
	class(calibrations) <- c(class(calibrations), "datelifeCalibrations")
	# instead of using structure()
	return(calibrations)
}

#' Extract secondary calibrations from a given \code{datelifeResult} object
#'
#' @description \code{extract_calibrations_dateliferesult} extracts divergence
#' times (i.e., secondary calibrations) for each taxon pair in a given
#' '\code{datelifeResult} object.
#'
#' @details The function calls summarize_datelife_result()] with
#' \code{summary_format = "phylo_all"} to go from a \code{datelifeResult} object
#' to a \code{phylo} or \code{multiPhylo} object that will be passed to
#' [extract_calibrations_phylo()].
#'
#' @param input A \code{datelifeResult} object.
#' @inheritParams get_all_calibrations
#' @inherit extract_calibrations_phylo return
#' @export
extract_calibrations_dateliferesult <- function(input = NULL,
																                each = FALSE) {

  phyloall <- suppressMessages(
              summarize_datelife_result(datelife_result = input,
                                        summary_format = "phylo_all"))

  return(extract_calibrations_phylo(input = phyloall,
                                    each = each))
}

#' Search and extract available secondary calibrations for a given character
#' vector of taxon names
#'
#' @description The function searches DateLife's local
#' database of phylogenetic trees with branch lengths proportional to time (aka,
#' chronograms) with [datelife_search()], and extracts divergence times
#' (i.e., secondary calibrations) from chronograms for each pair of given
#' taxon names with [extract_calibrations_phylo()].
#'
#' @details The function calls [datelife_search()]
#' with \code{summary_format = "phylo_all"} to get all chronograms in database
#' containing at least two taxa from \code{input}, and generates a \code{phylo}
#' or \code{multiPhylo} object object that will be passed to
#' [extract_calibrations_phylo()].
#'
#' @param input A character vector of taxon names.
#' @inheritParams get_all_calibrations
#' @inherit extract_calibrations_phylo return
#' @export
get_calibrations_vector <- function(input = NULL,
																    each = FALSE) {
  # TODO: is_datelife_search_input function or any type of input format checking
  # function to trap the case were input is a list
	phyloall <- datelife_search(input = input,
															summary_format = "phylo_all")

	return(extract_calibrations_phylo(input = phyloall,
                                    each = each))
}
#' Search and extract available secondary calibrations from a given
#' '\code{datelifeQuery} object
#'
#' @param input A \code{datelifeQuery} object.
#' @inheritParams get_all_calibrations
#' @inherit get_calibrations_vector description details
#' @inherit extract_calibrations_phylo return
#' @export
get_calibrations_datelifequery <- function(input = NULL,
																    			 each = FALSE) {
  if (suppressMessages(!is_datelife_query(input))) {
		stop("'input' is not a 'datelifeQuery' object.")
	}
	phyloall <- datelife_search(input = input,
															summary_format = "phylo_all")

	return(extract_calibrations_phylo(input = phyloall,
                                    each = each))
}
