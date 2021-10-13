#' Get secondary calibrations from a set of given taxon names
#'
#' \code{get_all_calibrations} performs a \code{\link[=datelife_search]{datelife_search}}
#' and gets divergence times (i.e., secondary calibrations) for each taxon name
#' pair given in \code{input}.
#'
#' @inheritParams datelife_search
# #'	\item{A \code{datelifeResult} object}{From \code{\link[=get_datelife_result]{get_datelife_result}}.}
#' @param each Boolean, default to \code{FALSE}: all calibrations are returned in
#' the same data frame. If \code{TRUE}, calibrations from each chronogram are returned
#' in separate data frames.
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
get_all_calibrations <- function(input = NULL,
																 each = FALSE) {
  if (inherits(input, "phylo") | inherits(input, "multiPhylo")) {
		if (inherits(input, "multiPhylo")) {
			tip_labels <- unname(unlist(lapply(input, "[", "tip.label")))
		} else {
			tip_labels <- input$tip.label
		}
		tip_labels <- unique(gsub(" ", "_", tip_labels))
    res <- get_calibrations_vector(input = tip_labels,
                                   each = each)
  }
  if (inherits(input, "character")) {
    res <- get_calibrations_vector(input = input,
                                   each = each)
  }
  if (inherits(input, "datelifeQuery")) {
    res <- get_calibrations_datelifequery(input = input,
                                   				 each = each)
  }
  return(res)
}
