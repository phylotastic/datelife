#' Get secondary calibrations from a set of given taxon names.
#'
#' \code{get_all_calibrations} gets divergence times (i.e., secondary
#' calibrations) for each taxon name pair given in \code{input}.
#'
#' @param input
#' \describe{
#'	\item{Taxon names}{As a character vector}
#'	\item{A tree with taxon names as tip labels}{As a \code{phylo} or \code{multiPhylo} object, OR as a newick character string}
#'	\item{A \code{datelifeResult} object}{From \code{\link[=get_datelife_result]{get_datelife_result}}
#'	}
#' @param each Boolean, default to \code{FALSE}: all calibrations are returned in
#' the same data frame. If \code{TRUE}, calibrations from each chronogram are returned
#' in separate data frames.
#' @return A data frame of secondary calibrations (or list of data frames, if \code{each = TRUE}),
#' for each pair of given taxon names. The attribute "chronograms" contains the
#' source chronograms from which the calibrations were obtained.
#' @export
get_all_calibrations <- function(input = NULL,
																 each = FALSE) {
  if(inherits(input, "phylo") | inherits(input, "multiPhylo")){
    res <- get_calibrations_phylo(input = input,
                                  each = each)
  }
  if(inherits(input, "character")){
    res <- get_calibrations_vector(input = input,
                                   each = each)
  }
  if(inherits(input, "datelifeResult")){
    res <- get_calibrations_vector(input = input,
                                   each = each)
  }
  return(res)
}
