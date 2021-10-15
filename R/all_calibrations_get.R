#' Extract secondary calibrations from a \code{phylo} or \code{multiPhylo} object with
#' branch lengths proportional to time.
#'
#' @description \code{extract_calibrations_phylo} extracts divergence times (i.e., secondary
#' calibrations) for each taxon pair in a given \code{phylo} or \code{multiPhylo}
#' object.
#'
#' @param phy a \code{phylo} or \code{multiPhylo} object with branch lengths proportional to time.
#' @inheritParams get_all_calibrations
#' @return A \code{data frame} (or list of data frames if \code{each = TRUE}) of secondary
#' calibrations available for each pair of taxon names given as tip labels in an
#' \code{input} \code{phylo}
#' or \code{multiPhylo} object. The attribute "chronograms" contains the \{input]
#' chronograms from which the secondary calibrations were extracted.
#'
#' @export
extract_calibrations_phylo <- function(phy = NULL,
																       each = FALSE) {
	if (inherits(phy, "multiPhylo")) {
		chronograms <- phy
		xx <- sapply(chronograms, "[", "edge.length")
		xx <- unname(sapply(xx, is.null))
	  if (all(xx)) {
	    message("trees in 'multiPhylo' input have no branch lengths. \n There are no calibrations to return.")
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
	if (inherits(phy, "phylo")) {
	  if (is.null(input$edge.length)) {
	    stop("'input' tree has no branch lengths.")
	  }
		chronograms <- list(phy)
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

#' Search and extract available secondary calibrations for a given character vector of taxon names.
#'
#' @description \code{get_calibrations_vector} searches DateLife's local database of phylogenetic
#' trees with branch lengths proportional to time (aka, chronograms) with
#' \code{\link[=datelife_search]{datelife_search}}, and extracts divergence times
#' (i.e., secondary calibrations) from chronograms for each pair of given taxon names.
#'
#' @details The function calls \code{\link[=datelife_search]{datelife_search}}
#' with \code{summary_format = "phylo_all"} to get all chronograms in database
#' containing at least two taxa from \code{input}, and generates a \code{phylo}
#' or \code{multiPhylo} object object that will be passed to
#' \code{\link[=get_calibrations_phylo]{get_calibrations_phylo}}.
#'
#' @param input A character vector of taxon names.
#' @inheritParams get_all_calibrations
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
get_calibrations_vector <- function(input = NULL,
																    each = FALSE) {
  # TODO: is_datelife_search_input function or any type of input format checking
  # function to trap the case were input is a list
	phyloall <- datelife_search(input = input,
															summary_format = "phylo_all")

	return(extract_calibrations_phylo(phy = phyloall,
                                    each = each))
}
#' Search and extract available secondary calibrations from a given \code{datelifeQuery} object.
#'
#' @description \code{get_calibrations_vector} searches DateLife's local database of phylogenetic
#' trees with branch lengths proportional to time (aka, chronograms) with
#' \code{\link[=datelife_search]{datelife_search}}, and extracts divergence times
#' (i.e., secondary calibrations) from chronograms for each pair of given taxon names.
#'
#' @details The function calls \code{\link[=datelife_search]{datelife_search}}
#' with \code{summary_format = "phylo_all"} to get all chronograms in database
#' containing at least two taxa from \code{input}, and generates a \code{phylo}
#' or \code{multiPhylo} object object that will be passed to
#' \code{\link[=get_calibrations_phylo]{get_calibrations_phylo}}.
#'
#' @param input A \code{datelifeQuery} object.
#' @inheritParams get_all_calibrations
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
get_calibrations_datelifequery <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
																    			 each = FALSE) {
  if (suppressMessages(!is_datelife_query(input))) {
		stop("'input' is not a 'datelifeQuery' object.")
	}
	phyloall <- datelife_search(input = input,
															summary_format = "phylo_all")

	return(extract_calibrations_phylo(phy = phyloall,
                                    each = each))
}
#' Extract secondary calibrations from a given \code{datelifeResult} object.
#'
#' @description \code{extract_calibrations_dateliferesult} extracts divergence times (i.e., secondary
#' calibrations) for each taxon pair in a given \code{datelifeResult} object.
#'
#' @details The function calls
#' \code{\link[=summarize_datelife_result]{summarize_datelife_result}} with
#' \code{summary_format = "phylo_all"} to go from a \code{datelifeResult} object
#' to a \code{phylo} or \code{multiPhylo} object that will be passed to
#' \code{\link[=get_calibrations_phylo]{get_calibrations_phylo}}.
#'
#' @param input A \code{datelifeResult} object.
#' @inheritParams get_all_calibrations
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
extract_calibrations_dateliferesult <- function(input = NULL,
																                each = FALSE) {

  phyloall <- suppressMessages(
              summarize_datelife_result(datelife_result = input,
                                        summary_format = "phylo_all"))

  return(extract_calibrations_phylo(phy = phyloall,
                                    each = each))
}
