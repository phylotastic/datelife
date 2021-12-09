#' Get secondary calibrations from a set of given taxon names
#'
#' @aliases datelife_calibrations
#' @description \code{get_all_calibrations} performs a [datelife_search()]
#' and gets divergence times (i.e., secondary calibrations) for each taxon name
#' pair given in \code{input}.
#'
#' @inheritParams datelife_search
#' @inheritParams extract_calibrations_phylo
#' @inherit extract_calibrations_phylo return
#' @export
get_all_calibrations <- function(input = NULL,
                                 each = FALSE) {
  #
  if (inherits(input, "datelifeQuery")) {
    res <- get_calibrations_datelifequery(
      datelife_query = input,
      each = each
    )
    return(res)
  }
  if (inherits(input, "phylo") | inherits(input, "multiPhylo")) {
    if (inherits(input, "multiPhylo")) {
      input <- unname(unlist(lapply(input, "[", "tip.label")))
    } else {
      input <- input$tip.label
    }
  }
  if (inherits(input, "character")) {
    input <- unique(gsub(" ", "_", input))
    res <- get_calibrations_vector(
      input = input,
      each = each
    )
  }
  return(res)
}
